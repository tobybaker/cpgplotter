"""
Input validation functions.

This module provides validation for genomic regions, BAM files, and other inputs.
"""

from pathlib import Path
import re
import pysam


def validate_region(region: str) -> tuple[str, int, int]:
    """
    Validate and parse genomic region string.

    Args:
        region: Region string in format "chr:start-end"

    Returns:
        Tuple of (chrom, start, end)

    Raises:
        ValueError: If region format is invalid

    Examples:
        >>> validate_region("chr7:1072064-1101499")
        ("chr7", 1072064, 1101499)
    """
    pattern = r"^([^:]+):(\d+)-(\d+)$"
    match = re.match(pattern, region)

    if not match:
        raise ValueError(
            f"Invalid region format: {region}. Expected 'chr:start-end'"
        )

    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))

    if start >= end:
        raise ValueError(f"Invalid region: start ({start}) must be < end ({end})")

    return chrom, start, end


def validate_bam_file(bam_path: Path) -> None:
    """
    Validate that BAM file exists and is indexed.

    Args:
        bam_path: Path to BAM file

    Raises:
        FileNotFoundError: If BAM or index file not found
        ValueError: If BAM file is invalid
    """
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    # Check for index (.bai or .csi)
    bai_path = Path(str(bam_path) + ".bai")
    csi_path = Path(str(bam_path) + ".csi")

    if not bai_path.exists() and not csi_path.exists():
        raise FileNotFoundError(
            f"BAM index not found for {bam_path}. Expected {bai_path} or {csi_path}"
        )

    # Try to open with pysam to validate
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # Try to get header to verify it's valid
            _ = bam.header
    except Exception as e:
        raise ValueError(f"Invalid BAM file {bam_path}: {e}")
