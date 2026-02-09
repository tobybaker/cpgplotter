"""
Side annotation axis rendering.

This module implements rendering of per-read annotation tracks adjacent to
the main heatmap. Qualitative annotations use discrete color palettes;
quantitative annotations use continuous colormaps.
"""

from typing import Literal

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch


def render_side_axis(
    ax: Axes,
    annotation_values: np.ndarray,
    annotation_type: Literal["qualitative", "quantitative"],
    palette: str | None = None,
    label: str = "",
) -> dict:
    """
    Render a side annotation axis.

    Dispatches to qualitative or quantitative rendering based on type.

    Args:
        ax: Matplotlib axis (narrow vertical strip, y-shared with heatmap).
        annotation_values: Array of annotation values (length = n_reads).
        annotation_type: "qualitative" or "quantitative".
        palette: Matplotlib colormap or palette name. Defaults to "tab10"
                 for qualitative, "viridis" for quantitative.
        label: Axis label text.

    Returns:
        Dict with legend info. For qualitative: {"handles": [...], "labels": [...]}.
        For quantitative: {"image": AxesImage}.
    """
    if annotation_type == "qualitative":
        default_palette = palette or "tab10"
        return _render_qualitative(ax, annotation_values, default_palette, label)
    else:
        default_palette = palette or "viridis"
        return _render_quantitative(ax, annotation_values, default_palette, label)


def _render_qualitative(
    ax: Axes,
    values: np.ndarray,
    palette: str,
    label: str,
) -> dict:
    """
    Render a qualitative (categorical) side axis.

    Each unique category gets a distinct color. Rendered as a (n_reads, 1) image.

    Args:
        ax: Matplotlib axis.
        values: Array of categorical values (strings or ints).
        palette: Matplotlib colormap name for category colors.
        label: Axis label.

    Returns:
        Dict with "handles" and "labels" for legend creation.
    """
    # Map categories to integer codes
    categories = sorted(set(values))
    cat_to_int = {cat: i for i, cat in enumerate(categories)}
    int_values = np.array([cat_to_int[v] for v in values], dtype=float)

    # Build discrete colormap
    n_cats = len(categories)
    base_cmap = plt.get_cmap(palette)
    colors = [base_cmap(i / max(n_cats - 1, 1)) for i in range(n_cats)]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(boundaries=np.arange(-0.5, n_cats), ncolors=n_cats)

    # Render as narrow image
    ax.imshow(
        int_values.reshape(-1, 1),
        aspect="auto",
        interpolation="none",
        cmap=cmap,
        norm=norm,
        origin="upper",
    )

    ax.set_xlabel(label, fontsize=8)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(axis="y", left=False, labelleft=False)

    # Build legend handles
    handles = [Patch(facecolor=colors[i], label=str(cat)) for i, cat in enumerate(categories)]
    return {"handles": handles, "labels": [str(c) for c in categories]}


def _render_quantitative(
    ax: Axes,
    values: np.ndarray,
    cmap_name: str,
    label: str,
) -> dict:
    """
    Render a quantitative (continuous) side axis.

    Rendered as a (n_reads, 1) image with a continuous colormap.

    Args:
        ax: Matplotlib axis.
        values: Array of numeric values.
        cmap_name: Matplotlib colormap name.
        label: Axis label.

    Returns:
        Dict with "image" for optional colorbar creation.
    """
    float_values = np.asarray(values, dtype=float)

    image = ax.imshow(
        float_values.reshape(-1, 1),
        aspect="auto",
        interpolation="none",
        cmap=cmap_name,
        origin="upper",
    )

    ax.set_xlabel(label, fontsize=8)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(axis="y", left=False, labelleft=False)

    return {"image": image}
