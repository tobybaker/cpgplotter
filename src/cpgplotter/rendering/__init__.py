"""
Visualization and rendering modules.

This package contains functions for rendering heatmaps, side axes, interval
overlays, and assembling multi-panel layouts.
"""

from cpgplotter.rendering.heatmap import render_heatmap, add_methylation_colorbar
from cpgplotter.rendering.side_axes import render_side_axis, render_qualitative_labels
from cpgplotter.rendering.overlays import render_interval_overlays
from cpgplotter.rendering.layout import (
    PanelAxes,
    PanelLayout,
    create_panel_layout,
    add_genomic_ticks,
)

__all__ = [
    "render_heatmap",
    "add_methylation_colorbar",
    "render_side_axis",
    "render_qualitative_labels",
    "render_interval_overlays",
    "PanelAxes",
    "PanelLayout",
    "create_panel_layout",
    "add_genomic_ticks",
]
