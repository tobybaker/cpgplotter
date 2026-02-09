"""
I/O utility functions.

This module provides helpers for loading sample sheets, parsing CLI arguments,
and handling file formats.
"""

from pathlib import Path
from typing import Optional
import polars as pl
from cpgplotter.core.config import SampleSpec


def load_sample_sheet(path: Path) -> list[SampleSpec]:
    """
    Load sample specifications from TSV file.

    Expected columns:
        - sample_name (required)
        - bam_path (required)
        - color (optional)
        - read_annotations (optional)
        - read_regions (optional)
        - Additional columns treated as metadata

    Args:
        path: Path to sample sheet TSV

    Returns:
        List of SampleSpec objects

    Example TSV format:
        sample_name    bam_path              color      read_annotations
        Tumor          /data/tumor.bam       #FC766A    tumor_annot.tsv
        Normal         /data/normal.bam      #184A45    .
    """
    df = pl.read_csv(path, separator="\t", comment_prefix="#")

    # Validate required columns
    required = {"sample_name", "bam_path"}
    if not required.issubset(set(df.columns)):
        raise ValueError(f"Sample sheet must have columns: {required}")

    samples = []
    for row in df.iter_rows(named=True):
        # Extract standard columns
        name = row["sample_name"]
        bam = Path(row["bam_path"])
        color = row.get("color")
        annot_path = row.get("read_annotations")
        regions_path = row.get("read_regions")

        # Handle "." as None
        if color == ".":
            color = None
        if annot_path in (".", None, ""):
            annot_path = None
        else:
            annot_path = Path(annot_path)
        if regions_path in (".", None, ""):
            regions_path = None
        else:
            regions_path = Path(regions_path)

        # Extract metadata (any other columns)
        metadata_keys = set(row.keys()) - {
            "sample_name",
            "bam_path",
            "color",
            "read_annotations",
            "read_regions",
        }
        metadata = {k: str(row[k]) for k in metadata_keys if row[k] is not None}

        samples.append(
            SampleSpec(
                name=name,
                bam=bam,
                color=color,
                read_annotations_path=annot_path,
                read_regions_path=regions_path,
                metadata=metadata,
            )
        )

    return samples


def parse_bam_spec(bam_spec: str) -> tuple[str, Path]:
    """
    Parse BAM specification from CLI argument.

    Format: "label:path" or just "path"
    If label is omitted, derives label from filename stem.

    Args:
        bam_spec: BAM specification string

    Returns:
        Tuple of (label, path)

    Examples:
        >>> parse_bam_spec("Tumor:/data/tumor.bam")
        ("Tumor", Path("/data/tumor.bam"))

        >>> parse_bam_spec("/data/tumor.bam")
        ("tumor", Path("/data/tumor.bam"))
    """
    if ":" in bam_spec:
        label, path_str = bam_spec.split(":", 1)
        return label, Path(path_str)
    else:
        path = Path(bam_spec)
        # Derive label from filename (without extensions)
        label = path.stem.replace(".sorted", "").replace(".bam", "")
        return label, path


def parse_figsize(figsize_str: str) -> tuple[float, float]:
    """
    Parse figure size string to tuple.

    Args:
        figsize_str: Figure size as "width,height" (e.g., "14,8")

    Returns:
        Tuple of (width, height) in inches

    Raises:
        ValueError: If format is invalid
    """
    try:
        parts = figsize_str.split(",")
        if len(parts) != 2:
            raise ValueError()
        width, height = float(parts[0]), float(parts[1])
        return width, height
    except (ValueError, AttributeError):
        raise ValueError(
            f"Invalid figsize format: {figsize_str}. Expected 'width,height'"
        )


def parse_annotation_spec(spec: str) -> dict[str, str]:
    """
    Parse annotation specification from CLI.

    Format: "label:source" where source is a BAM tag or "mapq"

    Args:
        spec: Annotation spec string

    Returns:
        Dictionary mapping label to source

    Example:
        >>> parse_annotation_spec("haplotype:HP")
        {"haplotype": "HP"}
    """
    if ":" not in spec:
        raise ValueError(f"Invalid annotation spec: {spec}. Expected 'label:source'")

    label, source = spec.split(":", 1)
    return {label: source}


def parse_side_axis_spec(spec: str) -> dict:
    """
    Parse side axis specification.

    Format: "column:type:palette" or "column:type" or "column"

    Args:
        spec: Side axis spec string

    Returns:
        Dictionary with keys: column, type (optional), palette (optional)

    Examples:
        >>> parse_side_axis_spec("haplotype:qualitative:Set2")
        {"column": "haplotype", "type": "qualitative", "palette": "Set2"}

        >>> parse_side_axis_spec("mapq")
        {"column": "mapq"}
    """
    parts = spec.split(":")
    result = {"column": parts[0]}

    if len(parts) >= 2:
        result["axis_type"] = parts[1]
    if len(parts) >= 3:
        result["palette"] = parts[2]

    return result
