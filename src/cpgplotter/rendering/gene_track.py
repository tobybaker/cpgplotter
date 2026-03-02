"""
Gene annotation track rendering.

This module renders gene models as a narrow, frameless axis showing
intron lines, exon boxes, strand arrows, and gene name labels — all
positioned in CpG coordinate space.
"""

from matplotlib.axes import Axes
from matplotlib.patches import FancyArrowPatch, Rectangle

from cpgplotter.core.coordinates import CpGIndex
from cpgplotter.core.gene_annotation import GeneModel, pack_gene_rows


def render_gene_track(
    ax: Axes,
    genes: list[GeneModel],
    cpg_index: CpGIndex,
    gene_color: str = "#2C3E50",
    label_fontsize: int = 7,
) -> None:
    """
    Render gene models on a frameless axis in CpG coordinate space.

    Each gene is drawn with a thin intron backbone line, thick exon
    rectangles, strand-indicating arrows, and a gene name label.

    Args:
        ax: Matplotlib axis to draw on
        genes: List of GeneModel objects (should already have representative
               transcripts selected)
        cpg_index: CpGIndex for genomic-to-CpG coordinate transformation
        gene_color: Color for gene model elements
        label_fontsize: Font size for gene name labels
    """
    rows = pack_gene_rows(genes)

    if not rows:
        _style_axis(ax)
        ax.set_ylim(0, 1)
        return

    row_height = 1.0
    row_spacing = 0.4
    total_height = len(rows) * (row_height + row_spacing)

    for row_idx, row in enumerate(rows):
        y_center = total_height - (row_idx + 0.5) * (row_height + row_spacing)

        for gene in row:
            if not gene.transcripts:
                continue
            tx = gene.transcripts[0]

            gene_start_cpg = cpg_index.genomic_to_cpg_continuous(gene.start)
            gene_end_cpg = cpg_index.genomic_to_cpg_continuous(gene.end)

            # Clamp to visible range
            vis_start = max(gene_start_cpg, -0.5)
            vis_end = min(gene_end_cpg, cpg_index.n_cpgs - 0.5)

            if vis_end <= vis_start:
                continue

            # Intron backbone line
            ax.plot(
                [vis_start, vis_end],
                [y_center, y_center],
                color=gene_color,
                linewidth=1,
                solid_capstyle="butt",
                zorder=1,
            )

            # Strand arrows along the intron line
            _draw_strand_arrows(
                ax, vis_start, vis_end, y_center,
                gene.strand, gene_color, row_height,
            )

            # Exon boxes
            exon_height = row_height * 0.6
            for exon_start, exon_end in tx.exons:
                ex_start_cpg = cpg_index.genomic_to_cpg_continuous(exon_start)
                ex_end_cpg = cpg_index.genomic_to_cpg_continuous(exon_end)

                # Clamp to visible range
                ex_start_cpg = max(ex_start_cpg, -0.5)
                ex_end_cpg = min(ex_end_cpg, cpg_index.n_cpgs - 0.5)

                if ex_end_cpg <= ex_start_cpg:
                    continue

                rect = Rectangle(
                    (ex_start_cpg, y_center - exon_height / 2),
                    ex_end_cpg - ex_start_cpg,
                    exon_height,
                    facecolor=gene_color,
                    edgecolor=gene_color,
                    linewidth=0.5,
                    zorder=2,
                )
                ax.add_patch(rect)

            # Gene name label
            label_x = (vis_start + vis_end) / 2
            label_y = y_center + row_height * 0.45
            ax.text(
                label_x,
                label_y,
                gene.gene_name,
                ha="center",
                va="bottom",
                fontsize=label_fontsize,
                fontstyle="italic",
                color=gene_color,
                zorder=3,
            )

    ax.set_xlim(-0.5, cpg_index.n_cpgs - 0.5)
    ax.set_ylim(0, total_height)
    _style_axis(ax)


def _draw_strand_arrows(
    ax: Axes,
    x_start: float,
    x_end: float,
    y: float,
    strand: str,
    color: str,
    row_height: float,
) -> None:
    """
    Draw small chevron arrows along the intron line to indicate strand.

    Args:
        ax: Matplotlib axis
        x_start: Start position in CpG space
        x_end: End position in CpG space
        y: Y-center of the gene row
        strand: "+" or "-"
        color: Arrow color
        row_height: Height of a gene row (for sizing arrows)
    """
    span = x_end - x_start
    if span < 2:
        return

    # Place arrows at regular intervals, not too dense
    n_arrows = max(1, min(int(span / 3), 8))
    arrow_size = row_height * 0.15

    for i in range(1, n_arrows + 1):
        x = x_start + i * span / (n_arrows + 1)
        if strand == "+":
            ax.plot(
                [x - arrow_size, x, x - arrow_size],
                [y + arrow_size, y, y - arrow_size],
                color=color,
                linewidth=0.7,
                solid_capstyle="round",
                zorder=1,
            )
        else:
            ax.plot(
                [x + arrow_size, x, x + arrow_size],
                [y + arrow_size, y, y - arrow_size],
                color=color,
                linewidth=0.7,
                solid_capstyle="round",
                zorder=1,
            )


def _style_axis(ax: Axes) -> None:
    """Remove all spines, ticks, and labels for a clean frameless look."""
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(
        axis="both",
        left=False, right=False, top=False, bottom=False,
        labelleft=False, labelright=False, labeltop=False, labelbottom=False,
    )
