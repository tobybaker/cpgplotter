"""
CpG Methylation Heatmap Visualization Tool

A Python package for visualizing read-level CpG methylation patterns from
long-read nanopore or PacBio sequencing data.
"""

__version__ = "0.1.0"

# High-level API exports
from cpgplotter.api import plot_methylation, plot_methylation_from_config

# Core data structures
from cpgplotter.core.coordinates import CpGIndex
from cpgplotter.core.annotations import ReadAnnotations, ReadIntervals
from cpgplotter.core.config import SampleSpec, PlotConfig

# Core functions
from cpgplotter.core.extraction import extract_methylation
from cpgplotter.core.gene_annotation import load_gene_annotations

# Rendering functions
from cpgplotter.rendering.heatmap import render_heatmap
from cpgplotter.rendering.layout import create_panel_layout, add_genomic_ticks
from cpgplotter.rendering.gene_track import render_gene_track

__all__ = [
    "__version__",
    "plot_methylation",
    "plot_methylation_from_config",
    "CpGIndex",
    "ReadAnnotations",
    "ReadIntervals",
    "SampleSpec",
    "PlotConfig",
    "extract_methylation",
    "load_gene_annotations",
    "render_heatmap",
    "render_gene_track",
    "create_panel_layout",
    "add_genomic_ticks",
]
