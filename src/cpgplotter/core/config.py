"""
Configuration and data model definitions using Pydantic.

This module defines all Pydantic models for configuration and validation,
including SampleSpec, PlotConfig, and related specifications.
"""

from pathlib import Path
from typing import Optional, Literal
from pydantic import BaseModel, Field


class SampleSpec(BaseModel):
    """
    Specification for a single BAM sample.

    Attributes:
        name: Human-readable label for the sample
        bam: Path to sorted, indexed BAM file
        color: Optional hex color for this sample (e.g., "#FC766A")
        read_annotations_path: Path to per-read scalar annotation TSV
        read_regions_path: Path to per-read interval BED file
        metadata: Additional sample-level metadata (e.g., group, treatment)
    """

    name: str
    bam: Path
    color: Optional[str] = None
    read_annotations_path: Optional[Path] = None
    read_regions_path: Optional[Path] = None
    metadata: dict[str, str] = Field(default_factory=dict)

    class Config:
        frozen = False


class SideAxisSpec(BaseModel):
    """
    Specification for a side annotation axis.

    Attributes:
        column: Name of annotation column to display
        axis_type: Type of annotation (qualitative or quantitative)
        palette: Matplotlib colormap or palette name
    """

    column: str
    axis_type: Literal["qualitative", "quantitative"]
    palette: Optional[str] = None

    class Config:
        frozen = False


class PlotConfig(BaseModel):
    """
    Top-level configuration for a cpgplotter figure.

    This model captures all parameters needed to generate a figure,
    including samples, region, sorting, rendering options, and output settings.

    Attributes:
        region: Genomic region string (chr:start-end)
        samples: List of sample specifications
        sort_by: List of annotation column names for read ordering
        side_axes: List of side axis specifications
        colormap: Matplotlib colormap for methylation probabilities
        panel_height_mode: How to size panels (uniform or proportional)
        max_reads: Maximum number of reads per panel
        min_mapq: Minimum mapping quality filter
        figsize: Optional figure size as (width, height) in inches
        output: Output file path
        output_format: Output format (png, svg, pdf)
        dpi: Output resolution
    """

    region: str
    samples: list[SampleSpec]
    sort_by: list[str] = Field(default_factory=list)
    side_axes: list[SideAxisSpec] = Field(default_factory=list)
    colormap: str = "RdBu_r"
    panel_height_mode: Literal["uniform", "proportional"] = "uniform"
    max_reads: Optional[int] = None
    min_mapq: int = 0
    min_cpg_coverage: int = 3
    min_cpgs_per_read: int = 5
    figsize: Optional[tuple[float, float]] = None
    output: Optional[Path] = None
    output_format: str = "png"
    dpi: int = 150

    class Config:
        frozen = False

    def to_yaml(self, path: Path) -> None:
        """
        Serialize configuration to YAML file.

        Args:
            path: Output YAML file path
        """
        import yaml

        with open(path, "w") as f:
            yaml.dump(self.model_dump(), f, default_flow_style=False)

    @classmethod
    def from_yaml(cls, path: Path) -> "PlotConfig":
        """
        Load configuration from YAML file.

        Args:
            path: Input YAML file path

        Returns:
            PlotConfig instance
        """
        import yaml

        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)
