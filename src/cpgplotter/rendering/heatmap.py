"""
Heatmap rendering functions.

This module implements the core methylation heatmap visualization using
matplotlib's imshow with masked arrays for NaN handling.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.image import AxesImage


def render_heatmap(
    matrix: np.ndarray,
    ax: Axes,
    colormap: str = "RdBu_r",
    vmin: float = 0.0,
    vmax: float = 1.0,
    missing_color: str = "#F0F0F0",
) -> AxesImage:
    """
    Render a methylation probability heatmap on a matplotlib axis.

    The matrix should already be sorted (caller handles sorting). NaN values
    are masked and rendered as the missing_color background.

    Args:
        matrix: Methylation probability matrix (n_reads x n_cpgs).
                Values in [0.0, 1.0], NaN for missing data.
        ax: Matplotlib axis to render on.
        colormap: Matplotlib colormap name.
        vmin: Minimum value for color mapping.
        vmax: Maximum value for color mapping.
        missing_color: Color for missing/uncovered CpGs.

    Returns:
        AxesImage object (for colorbar creation).
    """
    # Set background color for masked (NaN) cells
    ax.set_facecolor(missing_color)

    # Create masked array where NaN values are masked
    masked_matrix = np.ma.masked_invalid(matrix)

    # Render heatmap
    image = ax.imshow(
        masked_matrix,
        aspect="auto",
        interpolation="none",
        cmap=colormap,
        vmin=vmin,
        vmax=vmax,
        origin="upper",
    )

    return image


def add_methylation_colorbar(
    fig: plt.Figure,
    image: AxesImage,
    cbar_ax: Axes,
    label: str = "Methylation probability",
) -> None:
    """
    Add a methylation probability colorbar to a dedicated axis.

    Args:
        fig: Matplotlib figure.
        image: AxesImage from render_heatmap.
        cbar_ax: Dedicated axis for the colorbar.
        label: Colorbar label text.
    """
    cbar = fig.colorbar(image, cax=cbar_ax)
    cbar.set_label(label)
