"""
High-level Python API for cpgplotter.

This module provides the main `plot_methylation` function for programmatic use,
orchestrating the full pipeline from data extraction to rendering.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from cpgplotter.core.coordinates import CpGIndex
from cpgplotter.core.config import SampleSpec, PlotConfig, SideAxisSpec
from cpgplotter.core.extraction import extract_methylation, MethylationData
from cpgplotter.core.annotations import ReadAnnotations, ReadIntervals
from cpgplotter.core.gene_annotation import (
    load_gene_annotations,
    select_representative_transcripts,
)
from cpgplotter.processing.sorting import sort_reads
from cpgplotter.rendering.layout import create_panel_layout, add_genomic_ticks
from cpgplotter.rendering.heatmap import render_heatmap, add_methylation_colorbar
from cpgplotter.rendering.side_axes import render_side_axis, render_qualitative_labels
from cpgplotter.rendering.overlays import render_interval_overlays
from cpgplotter.rendering.gene_track import render_gene_track


def plot_methylation(
    region: str,
    bams: dict[str, str | Path] | str | Path,
    read_annotations: Optional[dict[str, str] | str | Path] = None,
    read_regions: Optional[str | Path] = None,
    sort_by: Optional[list[str]] = None,
    side_axes: Optional[dict[str, dict]] = None,
    colormap: str = "RdBu_r",
    panel_height_mode: str = "uniform",
    max_reads: Optional[int] = None,
    min_mapq: int = 0,
    min_cpg_coverage: int = 3,
    min_cpgs_per_read: int = 5,
    nan_weight: float = 0.5,
    gtf: Optional[str | Path] = None,
    gene_types: Optional[list[str]] = None,
    figsize: Optional[tuple[float, float]] = None,
    output: Optional[str | Path] = None,
    output_format: str = "png",
    dpi: int = 150,
) -> Figure:
    """
    Generate a CpG methylation heatmap visualization.

    This is the main high-level API function for creating methylation heatmaps.
    It handles all the processing steps from data extraction to rendering.

    Args:
        region: Genomic region string (e.g., "chr7:1072064-1101499")
        bams: BAM file specification. Can be:
            - dict: {"Label1": "path1.bam", "Label2": "path2.bam"}
            - str/Path: Single BAM file path
        read_annotations: Per-read scalar annotations. Can be:
            - dict: {"label": "tag_name"} for BAM tag extraction
            - str/Path: Path to TSV file
        read_regions: Path to per-read interval BED file
        sort_by: List of annotation column names for read ordering
        side_axes: Side axis specifications as dict:
            {"column_name": {"type": "qualitative", "palette": "Set2"}}
        colormap: Matplotlib colormap for methylation probabilities
        panel_height_mode: "uniform" or "proportional"
        max_reads: Maximum number of reads per panel
        min_mapq: Minimum mapping quality filter
        min_cpg_coverage: Minimum number of reads that must cover a CpG position
            for it to be included in the index. Default: 3.
        min_cpgs_per_read: Minimum number of CpGs a read must cover in the
            plotted region to be included. Default: 5.
        figsize: Figure size as (width, height) in inches
        output: Output file path (if None, figure is not saved)
        output_format: Output format ("png", "svg", "pdf")
        dpi: Output resolution

    Returns:
        Matplotlib Figure object
    """
    # 1. Normalize BAM inputs to list of SampleSpec
    samples = _normalize_bam_input(bams)

    # 2. Build shared CpG index from all BAMs
    bam_paths = [str(s.bam) for s in samples]
    cpg_index = CpGIndex(region=region, bam_path=bam_paths, min_cpg_coverage=min_cpg_coverage)

    # 2b. Load gene annotations if GTF provided
    gene_models = None
    if gtf is not None:
        gene_models = load_gene_annotations(
            gtf_path=gtf,
            chrom=cpg_index.chrom,
            start=cpg_index.start,
            end=cpg_index.end,
            gene_types=gene_types,
        )
        gene_models = select_representative_transcripts(gene_models)

    # 3. Extract and process data for each sample
    sample_data = []
    for sample in samples:
        # Extract methylation
        meth_data = extract_methylation(
            bam_path=str(sample.bam),
            cpg_index=cpg_index,
            min_mapq=min_mapq,
            max_reads=max_reads,
            min_cpgs_per_read=min_cpgs_per_read,
        )

        # Load annotations
        annot = _load_annotations(
            read_annotations=read_annotations,
            sample=sample,
            region=region,
        )

        # Load intervals
        intervals = _load_intervals(
            read_regions=read_regions,
            sample=sample,
            cpg_index=cpg_index,
        )

        # Sort reads
        sorted_indices, sorted_names = sort_reads(
            methylation_matrix=meth_data.matrix,
            read_names=meth_data.read_names,
            annotations=annot,
            sort_by=sort_by,
            nan_weight=nan_weight,
        )
        sorted_matrix = meth_data.matrix[sorted_indices]

        sample_data.append({
            "sample": sample,
            "matrix": sorted_matrix,
            "read_names": sorted_names,
            "annotations": annot,
            "intervals": intervals,
            "sorted_indices": sorted_indices,
        })

    # 4. Normalize side axes specs
    side_axis_specs = _normalize_side_axes(side_axes)
    n_side_axes = len(side_axis_specs)

    # 4b. Resolve side axis types for layout (qualitative axes get label columns)
    side_axis_types = _resolve_side_axis_types(
        side_axis_specs, sample_data[0]["annotations"] if sample_data else None
    )

    # 5. Create layout
    has_gene_track = gene_models is not None and len(gene_models) > 0
    panel_labels = [s["sample"].name for s in sample_data]
    panel_read_counts = [len(s["read_names"]) for s in sample_data]
    layout = create_panel_layout(
        n_sample_panels=len(samples),
        n_side_axes=n_side_axes,
        side_axis_types=side_axis_types,
        panel_height_mode=panel_height_mode,
        panel_read_counts=panel_read_counts,
        panel_labels=panel_labels,
        figsize=figsize,
        dpi=dpi,
        has_gene_track=has_gene_track,
    )

    # 5b. Render gene track if present
    if has_gene_track and layout.gene_track_ax is not None:
        render_gene_track(
            ax=layout.gene_track_ax,
            genes=gene_models,
            cpg_index=cpg_index,
        )

    # 6. Render each panel
    first_image = None
    all_legend_handles = []

    for panel_idx, (panel_axes, sdata) in enumerate(
        zip(layout.panels, sample_data)
    ):
        # Render heatmap
        image = render_heatmap(
            matrix=sdata["matrix"],
            ax=panel_axes.heatmap,
            colormap=colormap,
        )
        if first_image is None:
            first_image = image

        # Set panel label
        panel_axes.heatmap.set_ylabel(sdata["sample"].name)

        # Render side axes
        for sa_idx, sa_spec in enumerate(side_axis_specs):
            if sa_idx >= len(panel_axes.side_axes):
                break
            annot = sdata["annotations"]
            if annot is None:
                continue
            col_name = sa_spec["column"]
            if col_name not in annot.data.columns:
                continue

            # Get annotation values aligned with sorted read order
            annot_filtered = annot.filter_reads(set(sdata["read_names"]))
            name_to_val = dict(
                zip(
                    annot_filtered.data["read_name"].to_list(),
                    annot_filtered.data[col_name].to_list(),
                )
            )
            ordered_values = np.array(
                [name_to_val.get(n) for n in sdata["read_names"]]
            )

            ann_type = sa_spec.get(
                "type", annot.column_types.get(col_name, "qualitative")
            )
            palette = sa_spec.get("palette")

            render_info = render_side_axis(
                ax=panel_axes.side_axes[sa_idx],
                annotation_values=ordered_values,
                annotation_type=ann_type,
                palette=palette,
                label=col_name,
            )

            # Render inline labels for qualitative axes
            if (
                render_info.get("type") == "qualitative"
                and sa_idx < len(panel_axes.side_label_axes)
                and panel_axes.side_label_axes[sa_idx] is not None
            ):
                render_qualitative_labels(
                    label_ax=panel_axes.side_label_axes[sa_idx],
                    blocks=render_info["blocks"],
                )

        # Render interval overlays
        intervals = sdata["intervals"]
        if (
            intervals is not None
            and intervals.cpg_data is not None
            and len(intervals.cpg_data) > 0
        ):
            read_order = {
                name: idx for idx, name in enumerate(sdata["read_names"])
            }
            overlay_handles = render_interval_overlays(
                ax=panel_axes.heatmap,
                cpg_intervals=intervals.cpg_data,
                read_order=read_order,
            )
            if panel_idx == 0:
                all_legend_handles.extend(overlay_handles)

    # 7. Add shared colorbar
    if first_image is not None:
        add_methylation_colorbar(layout.fig, first_image, layout.cbar_ax)

    # 8. Add genomic ticks to bottom panel
    if layout.panels:
        bottom_ax = layout.panels[-1].heatmap
        add_genomic_ticks(bottom_ax, cpg_index)

    # 9. Add legend if we have handles (interval overlays only)
    if all_legend_handles:
        layout.fig.legend(
            handles=all_legend_handles,
            loc="upper right",
            fontsize=8,
            framealpha=0.8,
        )

    # 10. Save if output path provided
    if output is not None:
        layout.fig.savefig(
            output,
            format=output_format,
            dpi=dpi,
            bbox_inches="tight",
        )

    return layout.fig


def plot_methylation_from_config(config: PlotConfig | str | Path) -> Figure:
    """
    Generate a plot from a PlotConfig object or YAML file.

    Args:
        config: PlotConfig object or path to YAML config file

    Returns:
        Matplotlib Figure object
    """
    if isinstance(config, (str, Path)):
        config = PlotConfig.from_yaml(Path(config))

    # Convert config samples to bams dict
    bams = {s.name: s.bam for s in config.samples}

    # Convert side_axes specs
    sa = None
    if config.side_axes:
        sa = {}
        for spec in config.side_axes:
            entry = {"type": spec.axis_type}
            if spec.palette:
                entry["palette"] = spec.palette
            sa[spec.column] = entry

    return plot_methylation(
        region=config.region,
        bams=bams,
        sort_by=config.sort_by or None,
        side_axes=sa,
        colormap=config.colormap,
        panel_height_mode=config.panel_height_mode,
        max_reads=config.max_reads,
        min_mapq=config.min_mapq,
        min_cpg_coverage=config.min_cpg_coverage,
        min_cpgs_per_read=config.min_cpgs_per_read,
        nan_weight=config.nan_weight,
        gtf=config.gtf,
        gene_types=config.gene_types,
        figsize=config.figsize,
        output=config.output,
        output_format=config.output_format,
        dpi=config.dpi,
    )


def _normalize_bam_input(
    bams: dict[str, str | Path] | str | Path,
) -> list[SampleSpec]:
    """Convert various BAM input formats to a list of SampleSpec."""
    if isinstance(bams, (str, Path)):
        path = Path(bams)
        name = path.stem.replace(".sorted", "").replace(".bam", "")
        return [SampleSpec(name=name, bam=path)]
    elif isinstance(bams, dict):
        return [
            SampleSpec(name=name, bam=Path(path))
            for name, path in bams.items()
        ]
    else:
        raise ValueError(f"Unsupported bams type: {type(bams)}")


def _load_annotations(
    read_annotations: Optional[dict[str, str] | str | Path],
    sample: SampleSpec,
    region: str,
) -> Optional[ReadAnnotations]:
    """Load annotations from the appropriate source."""
    # Sample-specific annotation file takes precedence
    if sample.read_annotations_path is not None:
        return ReadAnnotations.from_tsv(sample.read_annotations_path)

    if read_annotations is None:
        return None

    if isinstance(read_annotations, (str, Path)):
        return ReadAnnotations.from_tsv(Path(read_annotations))
    elif isinstance(read_annotations, dict):
        return ReadAnnotations.from_bam_tags(
            bam_path=str(sample.bam),
            region=region,
            tags=read_annotations,
        )
    return None


def _load_intervals(
    read_regions: Optional[str | Path],
    sample: SampleSpec,
    cpg_index: CpGIndex,
) -> Optional[ReadIntervals]:
    """Load interval annotations and transform to CpG space."""
    path = None

    # Sample-specific file takes precedence
    if sample.read_regions_path is not None:
        path = sample.read_regions_path
    elif read_regions is not None:
        path = Path(read_regions)

    if path is None:
        return None

    intervals = ReadIntervals.from_tsv(path)
    intervals.transform_to_cpg_space(cpg_index)
    return intervals


def _normalize_side_axes(
    side_axes: Optional[dict[str, dict]],
) -> list[dict]:
    """Normalize side axes specification to a list of dicts."""
    if side_axes is None:
        return []

    specs = []
    for col_name, spec in side_axes.items():
        entry = {"column": col_name}
        if "type" in spec:
            entry["type"] = spec["type"]
        if "palette" in spec:
            entry["palette"] = spec["palette"]
        specs.append(entry)

    return specs


def _resolve_side_axis_types(
    side_axis_specs: list[dict],
    first_annotations: Optional[ReadAnnotations],
) -> list[str]:
    """Resolve side axis types from specs and annotation inference.

    Returns a list of "qualitative" or "quantitative" for each side axis.
    """
    types = []
    for spec in side_axis_specs:
        if "type" in spec:
            types.append(spec["type"])
        elif first_annotations is not None:
            col = spec["column"]
            types.append(
                first_annotations.column_types.get(col, "qualitative")
            )
        else:
            types.append("qualitative")
    return types
