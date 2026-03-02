"""
Gene annotation loading from GTF files.

This module handles parsing GENCODE/Ensembl GTF files to extract gene models
for rendering as annotation tracks. Supports both plain .gtf files (line-by-line
scan) and tabix-indexed .gtf.gz files (fast region queries via pysam).
"""

from dataclasses import dataclass, field
from pathlib import Path

import pysam


@dataclass
class TranscriptModel:
    """A single transcript with its exon structure."""

    transcript_id: str
    transcript_name: str
    start: int
    end: int
    exons: list[tuple[int, int]] = field(default_factory=list)
    is_canonical: bool = False


@dataclass
class GeneModel:
    """A gene with its transcripts."""

    gene_name: str
    gene_id: str
    gene_type: str
    chrom: str
    start: int
    end: int
    strand: str
    transcripts: list[TranscriptModel] = field(default_factory=list)


def _parse_gtf_attributes(attr_string: str) -> dict[str, str]:
    """
    Parse GTF attribute column into a dict.

    GTF attributes are semicolon-separated key-value pairs like:
      gene_id "ENSG00000223972.5"; gene_name "DDX11L1"; level 2;

    Args:
        attr_string: Raw attribute string from GTF column 9

    Returns:
        Dict mapping attribute keys to values (quotes stripped)
    """
    attrs: dict[str, str] = {}
    for item in attr_string.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        # Split on first space: key "value" or key value
        parts = item.split(" ", 1)
        if len(parts) == 2:
            key = parts[0]
            val = parts[1].strip().strip('"')
            # Handle repeated keys (e.g., multiple "tag" entries)
            if key == "tag":
                existing = attrs.get("_tags", "")
                attrs["_tags"] = f"{existing},{val}" if existing else val
            else:
                attrs[key] = val
    return attrs


def _parse_gtf_line(line: str) -> tuple[str, str, int, int, str, dict[str, str]] | None:
    """
    Parse a single GTF line into (chrom, feature, start, end, strand, attrs).

    Returns None for comment lines.
    """
    if line.startswith("#"):
        return None
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 9:
        return None
    chrom = fields[0]
    feature = fields[2]
    start = int(fields[3]) - 1  # GTF is 1-based, convert to 0-based
    end = int(fields[4])  # GTF end is inclusive, but we keep as-is for half-open
    strand = fields[6]
    attrs = _parse_gtf_attributes(fields[8])
    return chrom, feature, start, end, strand, attrs


def _overlaps(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Check if two intervals overlap."""
    return start1 < end2 and start2 < end1


def load_gene_annotations(
    gtf_path: str | Path,
    chrom: str,
    start: int,
    end: int,
    gene_types: list[str] | None = None,
) -> list[GeneModel]:
    """
    Load gene annotations from a GTF file for a given region.

    Supports both plain .gtf files (line-by-line scan with filtering) and
    tabix-indexed .gtf.gz files (fast region query via pysam).

    Args:
        gtf_path: Path to GTF or tabix-indexed GTF.gz file
        chrom: Chromosome name
        start: Region start (0-based)
        end: Region end (exclusive)
        gene_types: List of gene biotypes to include (e.g., ["protein_coding"]).
                   If None, defaults to ["protein_coding"].
                   Pass ["all"] or an empty list to include all types.

    Returns:
        List of GeneModel objects overlapping the region
    """
    if gene_types is None:
        gene_types = ["protein_coding"]
    filter_types = gene_types and "all" not in gene_types

    gtf_path = Path(gtf_path)

    # Collect parsed lines
    parsed_lines: list[tuple[str, int, int, str, dict[str, str]]] = []

    if gtf_path.suffix == ".gz" and Path(str(gtf_path) + ".tbi").exists():
        # Tabix-indexed: fast region query
        with pysam.TabixFile(str(gtf_path)) as tbx:
            for line in tbx.fetch(chrom, start, end):
                result = _parse_gtf_line(line)
                if result is None:
                    continue
                _, feature, f_start, f_end, strand, attrs = result
                if feature in ("gene", "transcript", "exon"):
                    parsed_lines.append((feature, f_start, f_end, strand, attrs))
    else:
        # Plain GTF: line-by-line scan with early filtering
        with open(gtf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                # Quick check: tab-split only first field for chromosome filter
                tab_pos = line.index("\t")
                line_chrom = line[:tab_pos]
                if line_chrom != chrom:
                    continue

                result = _parse_gtf_line(line)
                if result is None:
                    continue
                _, feature, f_start, f_end, strand, attrs = result
                if feature not in ("gene", "transcript", "exon"):
                    continue
                if not _overlaps(f_start, f_end, start, end):
                    continue
                parsed_lines.append((feature, f_start, f_end, strand, attrs))

    # Build gene models from parsed lines
    genes: dict[str, GeneModel] = {}  # gene_id -> GeneModel
    transcripts: dict[str, TranscriptModel] = {}  # transcript_id -> TranscriptModel
    transcript_to_gene: dict[str, str] = {}  # transcript_id -> gene_id

    for feature, f_start, f_end, strand, attrs in parsed_lines:
        gene_id = attrs.get("gene_id", "")
        gene_type = attrs.get("gene_type", "")

        if filter_types and gene_type not in gene_types:
            continue

        if feature == "gene":
            if gene_id not in genes:
                genes[gene_id] = GeneModel(
                    gene_name=attrs.get("gene_name", gene_id),
                    gene_id=gene_id,
                    gene_type=gene_type,
                    chrom=chrom,
                    start=f_start,
                    end=f_end,
                    strand=strand,
                )

        elif feature == "transcript":
            tx_id = attrs.get("transcript_id", "")
            if not tx_id:
                continue
            tags = attrs.get("_tags", "")
            is_canonical = "Ensembl_canonical" in tags
            tx = TranscriptModel(
                transcript_id=tx_id,
                transcript_name=attrs.get("transcript_name", tx_id),
                start=f_start,
                end=f_end,
                is_canonical=is_canonical,
            )
            transcripts[tx_id] = tx
            transcript_to_gene[tx_id] = gene_id

            # Ensure gene exists even without a separate gene line
            if gene_id not in genes:
                genes[gene_id] = GeneModel(
                    gene_name=attrs.get("gene_name", gene_id),
                    gene_id=gene_id,
                    gene_type=gene_type,
                    chrom=chrom,
                    start=f_start,
                    end=f_end,
                    strand=strand,
                )

        elif feature == "exon":
            tx_id = attrs.get("transcript_id", "")
            if tx_id in transcripts:
                transcripts[tx_id].exons.append((f_start, f_end))

    # Attach transcripts to genes
    for tx_id, tx in transcripts.items():
        gene_id = transcript_to_gene.get(tx_id)
        if gene_id and gene_id in genes:
            # Sort exons by start position
            tx.exons.sort(key=lambda e: e[0])
            genes[gene_id].transcripts.append(tx)

    return list(genes.values())


def select_representative_transcripts(genes: list[GeneModel]) -> list[GeneModel]:
    """
    Select one representative transcript per gene.

    Prefers the Ensembl_canonical transcript; falls back to the longest.

    Args:
        genes: List of GeneModel objects

    Returns:
        Same list with each gene's transcripts reduced to one representative
    """
    for gene in genes:
        if not gene.transcripts:
            continue
        # Prefer canonical transcript
        canonical = [t for t in gene.transcripts if t.is_canonical]
        if canonical:
            gene.transcripts = [canonical[0]]
        else:
            # Fall back to longest transcript
            longest = max(gene.transcripts, key=lambda t: t.end - t.start)
            gene.transcripts = [longest]
    return genes


def pack_gene_rows(genes: list[GeneModel]) -> list[list[GeneModel]]:
    """
    Pack genes into non-overlapping rows for stacked display.

    Uses a greedy interval scheduling approach: assign each gene to the
    first row where it doesn't overlap any existing gene.

    Args:
        genes: List of GeneModel objects, will be sorted by start position

    Returns:
        List of rows, where each row is a list of non-overlapping GeneModel objects
    """
    if not genes:
        return []

    sorted_genes = sorted(genes, key=lambda g: g.start)
    rows: list[list[GeneModel]] = []

    for gene in sorted_genes:
        placed = False
        for row in rows:
            # Check if gene overlaps any gene in this row
            last_in_row = row[-1]
            if gene.start >= last_in_row.end:
                row.append(gene)
                placed = True
                break
        if not placed:
            rows.append([gene])

    return rows
