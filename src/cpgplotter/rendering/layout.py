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
from typing import Literal, Optional

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


@dataclass
class PanelAxes:
    """Axes for a single sample panel."""

    heatmap: Axes
    side_axes: list[Axes] = field(default_factory=list)
    side_label_axes: list[Optional[Axes]] = field(default_factory=list)
    label: str = ""


@dataclass
class PanelLayout:
    """Complete layout for a multi-panel figure."""

    fig: Figure
    panels: list[PanelAxes] = field(default_factory=list)
    cbar_ax: Axes = None
    gene_track_ax: Optional[Axes] = None


def create_panel_layout(
    n_sample_panels: int,
    n_side_axes: int = 0,
    side_axis_types: list[str] | None = None,
    panel_height_mode: Literal["uniform", "proportional"] = "uniform",
    panel_read_counts: list[int] | None = None,
    panel_labels: list[str] | None = None,
    figsize: tuple[float, float] | None = None,
    dpi: int = 150,
    has_gene_track: bool = False,
    gene_track_height_ratio: float = 0.2,
) -> PanelLayout:
    """
    Create a multi-panel figure with shared x-axis and optional side axes.

    Args:
        n_sample_panels: Number of sample panels (one per BAM).
        n_side_axes: Number of side annotation axes per panel.
        side_axis_types: Type of each side axis ("qualitative" or "quantitative").
            Qualitative axes get a narrow label column to their left.
            If None, no label columns are created.
        panel_height_mode: How to size panels — "uniform" or "proportional".
        panel_read_counts: Read counts per panel (required if proportional).
        panel_labels: Labels for each panel (sample names).
        figsize: Figure size (width, height) in inches. Auto-calculated if None.
        dpi: Figure resolution.
        has_gene_track: Whether to include a gene annotation track above panels.
        gene_track_height_ratio: Height ratio for the gene track row relative
            to a uniform sample panel (default: 0.2).

    Returns:
        PanelLayout with figure, panel axes, colorbar axis, and optional
        gene track axis.
    """
    if panel_labels is None:
        panel_labels = [f"Panel {i}" for i in range(n_sample_panels)]

    # Calculate figure size if not provided
    if figsize is None:
        width = 12 + n_side_axes * 1.0
        height = max(4, 3 * n_sample_panels)
        if has_gene_track:
            height += 1.5
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

    # Add gene track row at the top if requested
    n_rows = n_sample_panels
    gene_track_row = None
    if has_gene_track:
        height_ratios = [gene_track_height_ratio] + height_ratios
        n_rows += 1
        gene_track_row = 0

    # Sub-grid for all rows (gene track + sample panels, stacked vertically)
    panel_gs = GridSpecFromSubplotSpec(
        nrows=n_rows,
        ncols=1,
        subplot_spec=outer_gs[0, 0],
        height_ratios=height_ratios,
        hspace=0.3,
    )

    # Colorbar axis spans the full height
    cbar_ax = fig.add_subplot(outer_gs[0, 1])

    # Build per-panel axes
    # Column structure depends on side axis types: qualitative axes get a
    # narrow label column to their left.
    side_ax_ratio = 1.5
    label_ratio = 0.8
    heatmap_ratio = 20

    # Build column layout: track which GridSpec column maps to what
    col_ratios = []
    # For each side axis index, store its GridSpec column for the strip
    strip_col_indices = []
    # For each side axis index, store its GridSpec column for the label (or None)
    label_col_indices = []

    _sa_types = side_axis_types or []

    for sa_idx in range(n_side_axes):
        if sa_idx < len(_sa_types) and _sa_types[sa_idx] == "qualitative":
            label_col_indices.append(len(col_ratios))
            col_ratios.append(label_ratio)
        else:
            label_col_indices.append(None)
        strip_col_indices.append(len(col_ratios))
        col_ratios.append(side_ax_ratio)

    heatmap_col_idx = len(col_ratios)
    col_ratios.append(heatmap_ratio)
    n_cols = len(col_ratios)

    panels = []
    prev_heatmap_ax = None
    gene_track_ax = None

    # Create gene track axis if requested
    if has_gene_track:
        gt_inner_gs = GridSpecFromSubplotSpec(
            nrows=1,
            ncols=n_cols,
            subplot_spec=panel_gs[0, 0],
            width_ratios=col_ratios,
            wspace=0.05,
        )
        gene_track_ax = fig.add_subplot(gt_inner_gs[0, heatmap_col_idx])
        prev_heatmap_ax = gene_track_ax

        # Hide all non-heatmap columns in the gene track row
        for col_idx in range(n_cols):
            if col_idx != heatmap_col_idx:
                empty_ax = fig.add_subplot(gt_inner_gs[0, col_idx])
                empty_ax.set_visible(False)

    # Sample panel offset (shifted by 1 if gene track present)
    panel_start_row = 1 if has_gene_track else 0

    for i in range(n_sample_panels):
        row_idx = panel_start_row + i
        inner_gs = GridSpecFromSubplotSpec(
            nrows=1,
            ncols=n_cols,
            subplot_spec=panel_gs[row_idx, 0],
            width_ratios=col_ratios,
            wspace=0.05,
        )

        # Create heatmap axis — share x with first panel's heatmap (or gene track)
        if prev_heatmap_ax is None:
            heatmap_ax = fig.add_subplot(inner_gs[0, heatmap_col_idx])
            prev_heatmap_ax = heatmap_ax
        else:
            heatmap_ax = fig.add_subplot(
                inner_gs[0, heatmap_col_idx], sharex=prev_heatmap_ax
            )

        # Create side axes and their label axes
        side_axes = []
        side_label_axes: list[Optional[Axes]] = []
        first_col_idx = strip_col_indices[0] if strip_col_indices else heatmap_col_idx

        for j in range(n_side_axes):
            strip_col = strip_col_indices[j]
            side_ax = fig.add_subplot(inner_gs[0, strip_col], sharey=heatmap_ax)
            side_ax.tick_params(
                axis="x", bottom=False, labelbottom=False
            )
            # Only show y-tick labels on the leftmost column
            if strip_col > first_col_idx:
                side_ax.tick_params(axis="y", left=False, labelleft=False)
            side_axes.append(side_ax)

            # Create label axis for qualitative side axes
            lbl_col = label_col_indices[j]
            if lbl_col is not None:
                lbl_ax = fig.add_subplot(inner_gs[0, lbl_col], sharey=heatmap_ax)
                lbl_ax.set_frame_on(False)
                lbl_ax.tick_params(
                    axis="both", left=False, labelleft=False,
                    bottom=False, labelbottom=False,
                )
                side_label_axes.append(lbl_ax)
            else:
                side_label_axes.append(None)

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
                side_label_axes=side_label_axes,
                label=panel_labels[i],
            )
        )

    return PanelLayout(
        fig=fig, panels=panels, cbar_ax=cbar_ax, gene_track_ax=gene_track_ax
    )


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
