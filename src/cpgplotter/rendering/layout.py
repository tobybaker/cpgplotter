"""
Multi-panel layout management.

This module implements functions for creating and arranging multi-panel figures
with shared axes and proper spacing. The layout uses a two-level GridSpec:
- Top level: one row per sample panel (+ colorbar column)
- Per panel: nested columns for side axes + heatmap

This design is extensible — additional track types (gene annotations, aggregate
plots) can be added as new rows in the top-level GridSpec without modifying
existing code.
"""

from dataclasses import dataclass, field
from typing import Literal

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


@dataclass
class PanelAxes:
    """Axes for a single sample panel."""

    heatmap: Axes
    side_axes: list[Axes] = field(default_factory=list)
    label: str = ""


@dataclass
class PanelLayout:
    """Complete layout for a multi-panel figure."""

    fig: Figure
    panels: list[PanelAxes] = field(default_factory=list)
    cbar_ax: Axes = None


def create_panel_layout(
    n_sample_panels: int,
    n_side_axes: int = 0,
    panel_height_mode: Literal["uniform", "proportional"] = "uniform",
    panel_read_counts: list[int] | None = None,
    panel_labels: list[str] | None = None,
    figsize: tuple[float, float] | None = None,
    dpi: int = 150,
) -> PanelLayout:
    """
    Create a multi-panel figure with shared x-axis and optional side axes.

    Args:
        n_sample_panels: Number of sample panels (one per BAM).
        n_side_axes: Number of side annotation axes per panel.
        panel_height_mode: How to size panels — "uniform" or "proportional".
        panel_read_counts: Read counts per panel (required if proportional).
        panel_labels: Labels for each panel (sample names).
        figsize: Figure size (width, height) in inches. Auto-calculated if None.
        dpi: Figure resolution.

    Returns:
        PanelLayout with figure, panel axes, and colorbar axis.
    """
    if panel_labels is None:
        panel_labels = [f"Panel {i}" for i in range(n_sample_panels)]

    # Calculate figure size if not provided
    if figsize is None:
        width = 12 + n_side_axes * 0.5
        height = max(4, 3 * n_sample_panels)
        figsize = (width, height)

    fig = plt.figure(figsize=figsize, dpi=dpi, layout="constrained")

    # Top-level grid: panels as rows, with a thin colorbar column on the right
    # Width ratios: [main_area, colorbar]
    outer_gs = GridSpec(
        nrows=1,
        ncols=2,
        figure=fig,
        width_ratios=[40, 1],
        wspace=0.05,
    )

    # Panel height ratios
    if panel_height_mode == "proportional" and panel_read_counts is not None:
        height_ratios = [max(c, 1) for c in panel_read_counts]
    else:
        height_ratios = [1] * n_sample_panels

    # Sub-grid for sample panels (stacked vertically)
    panel_gs = GridSpecFromSubplotSpec(
        nrows=n_sample_panels,
        ncols=1,
        subplot_spec=outer_gs[0, 0],
        height_ratios=height_ratios,
        hspace=0.3,
    )

    # Colorbar axis spans the full height
    cbar_ax = fig.add_subplot(outer_gs[0, 1])

    # Build per-panel axes
    # Column width ratios: [side_ax_1, side_ax_2, ..., heatmap]
    side_ax_ratio = 1
    heatmap_ratio = 20
    col_ratios = [side_ax_ratio] * n_side_axes + [heatmap_ratio]
    n_cols = n_side_axes + 1

    panels = []
    prev_heatmap_ax = None

    for i in range(n_sample_panels):
        inner_gs = GridSpecFromSubplotSpec(
            nrows=1,
            ncols=n_cols,
            subplot_spec=panel_gs[i, 0],
            width_ratios=col_ratios,
            wspace=0.05,
        )

        # Create heatmap axis — share x with first panel's heatmap
        if prev_heatmap_ax is None:
            heatmap_ax = fig.add_subplot(inner_gs[0, n_side_axes])
            prev_heatmap_ax = heatmap_ax
        else:
            heatmap_ax = fig.add_subplot(
                inner_gs[0, n_side_axes], sharex=prev_heatmap_ax
            )

        # Create side axes — share y with this panel's heatmap
        side_axes = []
        for j in range(n_side_axes):
            side_ax = fig.add_subplot(inner_gs[0, j], sharey=heatmap_ax)
            side_ax.tick_params(
                axis="x", bottom=False, labelbottom=False
            )
            # Only show y-tick labels on the leftmost axis
            if j > 0:
                side_ax.tick_params(axis="y", left=False, labelleft=False)
            side_axes.append(side_ax)

        # Hide x-tick labels on all but the bottom panel
        if i < n_sample_panels - 1:
            heatmap_ax.tick_params(axis="x", labelbottom=False)

        # Hide y-tick labels on heatmap if there are side axes
        if n_side_axes > 0:
            heatmap_ax.tick_params(axis="y", left=False, labelleft=False)

        panels.append(
            PanelAxes(
                heatmap=heatmap_ax,
                side_axes=side_axes,
                label=panel_labels[i],
            )
        )

    return PanelLayout(fig=fig, panels=panels, cbar_ax=cbar_ax)


def add_genomic_ticks(
    ax: Axes,
    cpg_index: object,
    n_ticks: int = 10,
) -> None:
    """
    Add genomic coordinate tick labels to axis.

    Args:
        ax: Matplotlib axis (typically bottom panel).
        cpg_index: CpGIndex object.
        n_ticks: Number of ticks to display.
    """
    tick_positions, tick_labels = cpg_index.get_tick_positions(n_ticks)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=45, ha="right")
    ax.set_xlabel(f"Genomic position ({cpg_index.chrom})")
