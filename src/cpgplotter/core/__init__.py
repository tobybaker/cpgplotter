"""
Core data structures and algorithms for cpgplotter.

This package contains the fundamental data structures (CpGIndex, ReadAnnotations,
ReadIntervals) and core algorithms for methylation extraction and coordinate
transformation.
"""

from cpgplotter.core.coordinates import CpGIndex
from cpgplotter.core.annotations import ReadAnnotations, ReadIntervals
from cpgplotter.core.config import SampleSpec, PlotConfig, SideAxisSpec
from cpgplotter.core.extraction import extract_methylation

__all__ = [
    "CpGIndex",
    "ReadAnnotations",
    "ReadIntervals",
    "SampleSpec",
    "PlotConfig",
    "SideAxisSpec",
    "extract_methylation",
]
