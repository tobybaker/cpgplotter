"""
Methylation data extraction from BAM files.

This module implements functions for extracting per-read, per-CpG methylation
probabilities from BAM files with MM/ML tags.
"""

from typing import Optional
import numpy as np
import pysam
from cpgplotter.core.coordinates import CpGIndex


class MethylationData:
    """
    Container for extracted methylation data.

    Attributes:
        read_names: Array of read names, shape (n_reads,)
        matrix: Methylation probability matrix, shape (n_reads, n_cpgs)
                Values are in [0.0, 1.0], NaN for missing/uncovered CpGs
        cpg_index: CpGIndex object defining the coordinate system
    """

    def __init__(
        self,
        read_names: np.ndarray,
        matrix: np.ndarray,
        cpg_index: CpGIndex,
    ):
        """
        Initialize MethylationData.

        Args:
            read_names: Array of read names
            matrix: Methylation probability matrix (n_reads x n_cpgs)
            cpg_index: CpGIndex for the region
        """
        if matrix.shape[0] != len(read_names):
            raise ValueError("Matrix rows must match number of read names")
        if matrix.shape[1] != cpg_index.n_cpgs:
            raise ValueError("Matrix columns must match number of CpGs")

        self.read_names = read_names
        self.matrix = matrix
        self.cpg_index = cpg_index

    @property
    def n_reads(self) -> int:
        """Number of reads in the dataset."""
        return len(self.read_names)

    @property
    def n_cpgs(self) -> int:
        """Number of CpGs in the region."""
        return self.cpg_index.n_cpgs

    def filter_by_mapq(self, bam_path: str, min_mapq: int) -> "MethylationData":
        """
        Filter reads by mapping quality.

        Args:
            bam_path: Path to BAM file (to read MAPQ values)
            min_mapq: Minimum mapping quality threshold

        Returns:
            New MethylationData with filtered reads
        """
        # TODO: Implement MAPQ filtering
        # 1. Open BAM and read MAPQ for each read
        # 2. Create boolean mask
        # 3. Filter read_names and matrix rows
        raise NotImplementedError("MAPQ filtering not yet implemented")

    def get_coverage(self) -> np.ndarray:
        """
        Calculate coverage (number of reads) per CpG.

        Returns:
            Array of coverage counts, shape (n_cpgs,)
        """
        return (~np.isnan(self.matrix)).sum(axis=0)

    def __repr__(self) -> str:
        return f"MethylationData(n_reads={self.n_reads}, n_cpgs={self.n_cpgs})"


def extract_methylation(
    bam_path: str,
    cpg_index: CpGIndex,
    min_mapq: int = 0,
    max_reads: Optional[int] = None,
    min_cpgs_per_read: int = 5,
) -> MethylationData:
    """
    Extract per-read methylation data from a BAM file.

    This function reads a BAM file, parses MM/ML tags for each read in the
    specified region, and constructs a methylation probability matrix aligned
    to the provided CpGIndex.

    Args:
        bam_path: Path to sorted, indexed BAM file
        cpg_index: CpGIndex defining the coordinate system
        min_mapq: Minimum mapping quality filter (default: 0)
        max_reads: Maximum number of reads to extract (default: all)
        min_cpgs_per_read: Minimum number of CpGs a read must cover in the
            plotted region to be included. Default: 5.

    Returns:
        MethylationData object containing read names and methylation matrix

    Implementation notes:
        - Only C+m (cytosine methylation) modifications are considered
        - MM tag encodes modification positions relative to reference
        - ML tag encodes modification probabilities (0-255 scale)
        - Missing CpGs (not covered by read) are represented as NaN
    """
    read_names_list = []
    matrix_rows = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(region=cpg_index.region):
            # Apply MAPQ filter
            if read.mapping_quality < min_mapq:
                continue

            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Parse methylation data from MM/ML tags
            positions, probabilities = parse_mm_ml_tags(read)

            # Initialize row with NaN (missing/uncovered CpGs)
            row = np.full(cpg_index.n_cpgs, np.nan, dtype=np.float64)

            # Map genomic positions to CpG indices and fill in probabilities
            for ref_pos, prob in zip(positions, probabilities):
                cpg_idx = cpg_index.genomic_to_cpg(ref_pos)
                if cpg_idx is not None:
                    row[cpg_idx] = prob

            # Skip reads with fewer than min_cpgs_per_read CpGs in the region
            if np.sum(~np.isnan(row)) < min_cpgs_per_read:
                continue

            read_names_list.append(read.query_name)
            matrix_rows.append(row)

            # Apply max_reads limit if specified
            if max_reads is not None and len(read_names_list) >= max_reads:
                break

    # Convert lists to numpy arrays
    if len(read_names_list) == 0:
        # No reads found - return empty data
        read_names = np.array([], dtype=object)
        matrix = np.zeros((0, cpg_index.n_cpgs), dtype=np.float64)
    else:
        read_names = np.array(read_names_list, dtype=object)
        matrix = np.array(matrix_rows, dtype=np.float64)

    return MethylationData(read_names, matrix, cpg_index)


def parse_mm_ml_tags(read: pysam.AlignedSegment) -> tuple[list[int], list[float]]:
    """
    Parse MM and ML tags from a BAM read.

    Args:
        read: pysam AlignedSegment

    Returns:
        Tuple of (positions, probabilities)
        - positions: List of genomic positions with C+m modifications
        - probabilities: List of methylation probabilities [0.0, 1.0]

    Raises:
        ValueError: If MM/ML tags are missing or malformed
    """
    if not hasattr(read, 'modified_bases') or read.modified_bases is None:
        return [], []

    positions = []
    probabilities = []

    # Look for C+m modifications (cytosine methylation)
    # The key format is (base, strand, modification_type)
    # For PacBio, we need to check both strands: ('C', 0, 'm') and ('C', 1, 'm')
    for mod_key, mod_positions in read.modified_bases.items():
        base, strand, mod_type = mod_key

        # Only process cytosine methylation
        if base == 'C' and mod_type == 'm':
            # Get reference positions for this read (maps query pos -> ref pos)
            ref_positions = read.get_reference_positions(full_length=True)

            for query_pos, prob_uint8 in mod_positions:
                # Convert query position to reference position
                if query_pos < len(ref_positions) and ref_positions[query_pos] is not None:
                    ref_pos = ref_positions[query_pos]
                    # Normalize reverse-strand CpG positions to the C coordinate
                    if strand == 1:
                        ref_pos -= 1

                    # Convert probability from 0-255 scale to 0.0-1.0
                    prob = prob_uint8 / 255.0

                    positions.append(ref_pos)
                    probabilities.append(prob)

    return positions, probabilities
