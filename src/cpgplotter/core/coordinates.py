"""
CpG coordinate system management.

This module implements the CpGIndex class, which maintains the canonical coordinate
system for a genomic region. The CpGIndex maps between genomic positions and
CpG indices, ensuring consistent coordinate transformation across all panels.
"""

from collections import Counter
from typing import Optional
import numpy as np
import polars as pl
import pysam


class CpGIndex:
    """
    Canonical coordinate system for a genomic region.

    The CpGIndex represents the ordered sequence of CpG dinucleotides within
    a region, providing bidirectional mapping between genomic positions and
    CpG indices. This index is constructed once and shared across all panels
    to guarantee coordinate consistency.

    Attributes:
        region: Original region string (e.g., "chr7:1072064-1101499")
        chrom: Chromosome name
        start: Region start position (0-based)
        end: Region end position (exclusive)
        positions: Genomic positions of each CpG, shape (N,)
        n_cpgs: Number of CpGs in the region
    """

    def __init__(
        self,
        region: str,
        bam_path: Optional[str | list[str]] = None,
        positions: Optional[np.ndarray] = None,
        min_cpg_coverage: int = 3,
    ):
        """
        Initialize CpG coordinate system for a region.

        Args:
            region: Genomic region string in format "chr:start-end"
            bam_path: Path(s) to BAM file(s) for extracting CpG positions from MM tags.
                     Can be a single path string or a list of paths. When multiple BAMs
                     are provided, CpG positions are collected from all BAMs to create
                     a shared coordinate system.
            positions: Pre-computed CpG positions array (alternative to bam_path)
            min_cpg_coverage: Minimum number of reads that must cover a CpG position
                            for it to be included in the index. Filters out spurious
                            CpG sites caused by sequencing errors. Default: 3.
                            Set to 1 to disable filtering.

        Either bam_path or positions must be provided.
        """
        self.region = region
        self.chrom, coords = region.split(":")
        self.start, self.end = map(int, coords.split("-"))
        self.min_cpg_coverage = min_cpg_coverage

        if positions is not None:
            self.positions = positions
        elif bam_path is not None:
            # Normalize to list for uniform handling
            bam_paths: list[str] = [bam_path] if isinstance(bam_path, str) else bam_path
            self.positions = self._extract_cpg_positions(bam_paths)
        else:
            raise ValueError("Either bam_path or positions must be provided")

        self.n_cpgs = len(self.positions)

    def _extract_cpg_positions(self, bam_paths: list[str]) -> np.ndarray:
        """
        Extract CpG positions from one or more BAM files.

        Reads the specified region from all provided BAM files and identifies
        CpG positions that have modification calls (MM/ML tags) in at least
        min_cpg_coverage reads. This filters out spurious CpG sites caused
        by sequencing errors.

        Args:
            bam_paths: List of paths to sorted, indexed BAM files

        Returns:
            Sorted array of CpG genomic positions meeting coverage threshold
        """
        cpg_coverage: Counter[int] = Counter()

        for bam_path in bam_paths:
            with pysam.AlignmentFile(bam_path, "rb") as bam:
                for read in bam.fetch(region=self.region):
                    # Skip unmapped reads
                    if read.is_unmapped:
                        continue

                    # Check for modified bases
                    if not hasattr(read, 'modified_bases') or read.modified_bases is None:
                        continue

                    # Collect unique CpG positions for this read
                    read_cpg_positions: set[int] = set()

                    # Look for C+m modifications (cytosine methylation)
                    for mod_key, mod_positions in read.modified_bases.items():
                        base, strand, mod_type = mod_key

                        # Only process cytosine methylation
                        if base == 'C' and mod_type == 'm':
                            # Get reference positions for this read
                            ref_positions = read.get_reference_positions(full_length=True)

                            for query_pos, _ in mod_positions:
                                # Convert query position to reference position
                                if query_pos < len(ref_positions) and ref_positions[query_pos] is not None:
                                    ref_pos = ref_positions[query_pos]
                                    # Normalize reverse-strand CpG positions to the C coordinate
                                    if strand == 1:
                                        ref_pos -= 1

                                    # Only include positions within our region
                                    if self.start <= ref_pos < self.end:
                                        read_cpg_positions.add(ref_pos)

                    # Count each position once per read
                    for pos in read_cpg_positions:
                        cpg_coverage[pos] += 1

        # Filter by minimum coverage threshold
        passing_positions = [
            pos for pos, count in cpg_coverage.items()
            if count >= self.min_cpg_coverage
        ]

        if len(passing_positions) == 0:
            raise ValueError(
                f"No CpG positions found in region {self.region} "
                f"with coverage >= {self.min_cpg_coverage}"
            )

        cpg_positions = np.array(sorted(passing_positions), dtype=np.int64)
        return cpg_positions

    def genomic_to_cpg(self, genomic_pos: int) -> Optional[int]:
        """
        Convert genomic position to CpG index.

        Args:
            genomic_pos: Genomic coordinate

        Returns:
            CpG index (0-based) if position corresponds to a CpG, None otherwise
        """
        idx = np.searchsorted(self.positions, genomic_pos)
        if idx < len(self.positions) and self.positions[idx] == genomic_pos:
            return int(idx)
        return None

    def cpg_to_genomic(self, cpg_idx: int) -> int:
        """
        Convert CpG index to genomic position.

        Args:
            cpg_idx: CpG index (0-based)

        Returns:
            Genomic position of the CpG

        Raises:
            IndexError: If cpg_idx is out of range
        """
        return int(self.positions[cpg_idx])

    def transform_intervals(self, intervals_df: pl.DataFrame) -> pl.DataFrame:
        """
        Transform interval coordinates from genomic to CpG space.

        Takes a DataFrame with genomic interval coordinates (start, end) and
        adds corresponding CpG space coordinates (cpg_start, cpg_end).

        For interval [g_start, g_end) in genomic coordinates, finds all CpG
        indices i where g_start <= positions[i] < g_end. The CpG-space interval
        is [min(i), max(i)] (inclusive). Intervals with no CpGs are dropped.

        Args:
            intervals_df: DataFrame with columns: read_name, start, end, [label], [score]
                         start/end are genomic coordinates (0-based, half-open)

        Returns:
            DataFrame with added columns: cpg_start, cpg_end
            Intervals with no CpGs are filtered out
        """
        starts = intervals_df["start"].to_numpy()
        ends = intervals_df["end"].to_numpy()

        # Find first CpG index >= g_start
        cpg_starts = np.searchsorted(self.positions, starts, side="left")
        # Find last CpG index < g_end
        cpg_ends = np.searchsorted(self.positions, ends, side="left") - 1

        # Valid intervals have at least one CpG: cpg_start <= cpg_end
        # and indices must be within bounds
        valid = (cpg_starts <= cpg_ends) & (cpg_starts < self.n_cpgs) & (cpg_ends >= 0)

        result = intervals_df.with_columns(
            pl.Series("cpg_start", cpg_starts.astype(np.int64)),
            pl.Series("cpg_end", cpg_ends.astype(np.int64)),
            pl.Series("_valid", valid),
        ).filter(pl.col("_valid")).drop("_valid")

        return result

    def get_tick_positions(
        self, n_ticks: int = 10
    ) -> tuple[np.ndarray, list[str]]:
        """
        Generate tick positions and labels for genomic coordinate axis.

        Creates evenly-spaced tick marks in CpG space with corresponding
        genomic coordinate labels.

        Args:
            n_ticks: Approximate number of ticks to generate

        Returns:
            Tuple of (tick_positions, tick_labels)
            - tick_positions: Array of CpG indices for tick placement
            - tick_labels: List of formatted genomic coordinate strings
        """
        # Generate evenly spaced indices
        indices = np.linspace(0, self.n_cpgs - 1, n_ticks, dtype=int)

        # Get genomic positions for these indices
        genomic_positions = self.positions[indices]

        # Format labels (e.g., "1,072,064")
        labels = [f"{pos:,}" for pos in genomic_positions]

        return indices, labels

    def __repr__(self) -> str:
        return f"CpGIndex(region='{self.region}', n_cpgs={self.n_cpgs})"
