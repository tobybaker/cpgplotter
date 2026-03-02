"""
Interval overlay rendering.

This module implements rendering of per-read interval annotations as
semi-transparent rectangles overlaid on the heatmap.
"""

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle, Patch


def render_interval_overlays(
    ax: Axes,
    cpg_intervals: pl.DataFrame,
    read_order: dict[str, int],
    alpha: float = 0.8,
    label_colors: dict[str, str] | None = None,
) -> list[Patch]:
    """
    Render per-read interval overlays on a heatmap.

    Intervals are drawn as semi-transparent rectangles in CpG coordinate space.
    Different label values get different colors; score values modulate opacity.

    Args:
        ax: Matplotlib axis (the heatmap axis).
        cpg_intervals: DataFrame with columns: read_name, cpg_start, cpg_end,
                       and optionally label, score. Coordinates are in CpG space.
        read_order: Mapping of read_name -> y-position in plot.
        alpha: Base transparency for rectangles.
        label_colors: Optional mapping of label -> color. If None, colors are
                      assigned from a default palette.

    Returns:
        List of legend Patch handles for the interval labels.
    """
    if cpg_intervals is None or len(cpg_intervals) == 0:
        return []

    has_label = "label" in cpg_intervals.columns
    has_score = "score" in cpg_intervals.columns

    # Assign colors to labels
    if has_label:
        labels = cpg_intervals["label"].unique().sort().to_list()
    else:
        labels = ["interval"]

    if label_colors is None:
        # Custom bright palette: Lime Green, Yellow, Pink, Orange, Cyan, Bright White
        bright_colors = ["#32CD32", "#FFFF00", "#FF69B4", "#FFA500", "#00FFFF", "#FFFFFF"]
        
        label_colors = {
            lab: bright_colors[i % len(bright_colors)] for i, lab in enumerate(labels)
        }

    # Draw rectangles
    for row in cpg_intervals.iter_rows(named=True):
        read_name = row["read_name"]
        if read_name not in read_order:
            continue

        y_pos = read_order[read_name]
        cpg_start = row["cpg_start"]
        cpg_end = row["cpg_end"]

        label = row.get("label", "interval") if has_label else "interval"
        color = label_colors.get(label, "#888888")

        # Score modulates opacity
        rect_alpha = alpha

        width = cpg_end - cpg_start + 1
        rect = Rectangle(
            xy=(cpg_start - 0.5, y_pos - 0.5),
            width=width,
            height=1,
            facecolor=color,
            edgecolor=color,
            alpha=rect_alpha,
            linewidth=0.5,
        )
        ax.add_patch(rect)

    # Build legend handles
    handles = [
        Patch(facecolor=label_colors[lab], alpha=alpha, label=str(lab))
        for lab in labels
    ]
    return handles
