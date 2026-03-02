"""
Command-line interface for cpgplotter.

This module implements the Click-based CLI with the main 'plot' command and
associated options for visualizing CpG methylation patterns.
"""

import click
from pathlib import Path
from typing import Optional

from cpgplotter.utils.io import (
    parse_bam_spec,
    parse_figsize,
    parse_annotation_spec,
    parse_side_axis_spec,
    load_sample_sheet,
)
from cpgplotter.utils.validation import validate_region, validate_bam_file
from cpgplotter.core.config import SampleSpec, PlotConfig, SideAxisSpec


@click.group()
@click.version_option(version="0.1.0")
def main():
    """
    CpG Methylation Heatmap Visualisation Tool

    Visualize read-level CpG methylation patterns from long-read sequencing data.
    """
    pass


@main.command()
@click.option(
    "--region",
    type=str,
    help="Genomic region in format chr:start-end (e.g., chr7:1072064-1101499)",
)
@click.option(
    "--bam",
    "bams",
    type=str,
    multiple=True,
    help="BAM file as 'label:path' or just 'path'. Can be specified multiple times.",
)
@click.option(
    "--samples",
    type=click.Path(exists=True, path_type=Path),
    help="Sample sheet TSV file (mutually exclusive with --bam)",
)
@click.option(
    "--read-annotations",
    "read_annotations",
    type=str,
    multiple=True,
    help="Per-read scalar annotation as 'label:source' (BAM tag) or path to TSV",
)
@click.option(
    "--read-regions",
    "read_regions",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="Per-read interval BED file",
)
@click.option(
    "--sort-by",
    type=str,
    help="Comma-separated list of annotation columns for read ordering",
)
@click.option(
    "--side-axes",
    "side_axes",
    type=str,
    multiple=True,
    help="Side axis spec: 'column:type:palette' (e.g., haplotype:qualitative:Set2)",
)
@click.option(
    "--colormap",
    type=str,
    default="RdBu_r",
    help="Matplotlib colormap for methylation probabilities",
)
@click.option(
    "--panel-height-mode",
    type=click.Choice(["uniform", "proportional"]),
    default="uniform",
    help="Panel height mode: uniform or proportional to read count",
)
@click.option(
    "--min-mapq",
    type=int,
    default=0,
    help="Minimum mapping quality filter",
)
@click.option(
    "--min-cpg-coverage",
    type=int,
    default=3,
    help="Minimum number of reads covering a CpG site for it to be included (filters sequencing errors)",
)
@click.option(
    "--min-cpgs-per-read",
    type=int,
    default=5,
    help="Minimum number of CpGs a read must cover in the region to be included",
)
@click.option(
    "--nan-weight",
    type=float,
    default=0.5,
    help="Max weight for NaN coverage pattern in clustering distance (0=pure Hellinger, 1=coverage can dominate). Default: 0.5",
)
@click.option(
    "--max-reads",
    type=int,
    help="Maximum number of reads per panel",
)
@click.option(
    "--figsize",
    type=str,
    help="Figure size as 'width,height' in inches",
)
@click.option(
    "--dpi",
    type=int,
    default=150,
    help="Output resolution (dots per inch)",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=Path),
    help="Output file path",
)
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["png", "svg", "pdf"]),
    default="png",
    help="Output format",
)
@click.option(
    "--dump-config",
    type=click.Path(path_type=Path),
    help="Dump resolved configuration to YAML file",
)
@click.option(
    "--config",
    type=click.Path(exists=True, path_type=Path),
    help="Load configuration from YAML file (overrides other options)",
)
@click.option(
    "--gtf",
    type=click.Path(exists=True, path_type=Path),
    help="GTF annotation file (.gtf or tabix-indexed .gtf.gz) for gene track above heatmap",
)
@click.option(
    "--gene-types",
    type=str,
    default="protein_coding",
    help="Comma-separated gene biotypes to display (default: protein_coding). Use 'all' for all types.",
)
def plot(
    region: Optional[str],
    bams: tuple[str, ...],
    samples: Optional[Path],
    read_annotations: tuple[str, ...],
    read_regions: tuple[Path, ...],
    sort_by: Optional[str],
    side_axes: tuple[str, ...],
    colormap: str,
    panel_height_mode: str,
    min_mapq: int,
    min_cpg_coverage: int,
    min_cpgs_per_read: int,
    nan_weight: float,
    max_reads: Optional[int],
    figsize: Optional[str],
    dpi: int,
    output: Optional[Path],
    output_format: str,
    dump_config: Optional[Path],
    config: Optional[Path],
    gtf: Optional[Path],
    gene_types: str,
):
    """
    Generate a CpG methylation heatmap visualization.

    This command creates a heatmap showing per-read CpG methylation patterns
    across a specified genomic region.

    Examples:

        # Single BAM file
        cpgplotter plot --bam tumor.bam --region chr7:1072064-1101499

        # Multiple BAMs with annotations
        cpgplotter plot \\
            --bam Tumor:tumor.bam \\
            --bam Normal:normal.bam \\
            --region chr7:1072064-1101499 \\
            --read-annotations haplotype:HP \\
            --sort-by haplotype \\
            --side-axes haplotype:qualitative:Set2
    """
    from cpgplotter.api import plot_methylation, plot_methylation_from_config

    # Handle --config replay mode
    if config is not None:
        click.echo(f"Loading configuration from {config}")
        fig = plot_methylation_from_config(config)
        click.echo("Plot generated successfully.")
        return

    # Validate required inputs
    if not bams and not samples:
        raise click.UsageError("Must provide either --bam, --samples, or --config")

    if bams and samples:
        raise click.UsageError("--bam and --samples are mutually exclusive")

    if region is None:
        raise click.UsageError("--region is required unless using --config")

    # Validate region format
    validate_region(region)

    # Build BAMs dict or sample specs
    if samples:
        sample_specs = load_sample_sheet(samples)
        bams_arg = {s.name: s.bam for s in sample_specs}
    else:
        bams_arg = {}
        for bam_spec in bams:
            label, path = parse_bam_spec(bam_spec)
            validate_bam_file(path)
            bams_arg[label] = path

    # Parse read annotations
    annot_arg = None
    if read_annotations:
        # Check if first annotation looks like a file path
        first = read_annotations[0]
        if Path(first).exists() and not ":" in first:
            annot_arg = first
        else:
            # Parse as label:source BAM tag specs
            tags = {}
            for spec in read_annotations:
                if Path(spec).exists():
                    # It's a TSV file path — use it directly
                    annot_arg = spec
                    break
                tags.update(parse_annotation_spec(spec))
            if tags:
                annot_arg = tags

    # Parse read regions (use first one as global)
    regions_arg = None
    if read_regions:
        regions_arg = read_regions[0]

    # Parse sort-by
    sort_by_list = None
    if sort_by:
        sort_by_list = [s.strip() for s in sort_by.split(",")]

    # Parse side axes
    side_axes_dict = None
    if side_axes:
        side_axes_dict = {}
        for spec_str in side_axes:
            parsed = parse_side_axis_spec(spec_str)
            col = parsed.pop("column")
            side_axes_dict[col] = parsed

    # Parse figsize
    figsize_tuple = None
    if figsize:
        figsize_tuple = parse_figsize(figsize)

    # Auto-generate output path if not specified
    if output is None:
        region_clean = region.replace(":", "_").replace("-", "_")
        output = Path(f"cpgplotter_{region_clean}.{output_format}")

    # Parse gene types
    if gene_types == "all":
        gene_types_list = ["all"]
    else:
        gene_types_list = [t.strip() for t in gene_types.split(",")]

    click.echo(f"Plotting region: {region}")
    click.echo(f"Samples: {len(bams_arg)}")
    if gtf:
        click.echo(f"Gene annotation: {gtf}")

    # Generate plot
    fig = plot_methylation(
        region=region,
        bams=bams_arg,
        read_annotations=annot_arg,
        read_regions=regions_arg,
        sort_by=sort_by_list,
        side_axes=side_axes_dict,
        colormap=colormap,
        panel_height_mode=panel_height_mode,
        max_reads=max_reads,
        min_mapq=min_mapq,
        min_cpg_coverage=min_cpg_coverage,
        min_cpgs_per_read=min_cpgs_per_read,
        nan_weight=nan_weight,
        gtf=gtf,
        gene_types=gene_types_list,
        figsize=figsize_tuple,
        output=output,
        output_format=output_format,
        dpi=dpi,
    )

    click.echo(f"Plot saved to: {output}")

    # Dump config if requested
    if dump_config is not None:
        sample_list = [
            SampleSpec(name=name, bam=Path(path))
            for name, path in bams_arg.items()
        ]
        sa_specs = []
        if side_axes_dict:
            for col, spec in side_axes_dict.items():
                sa_specs.append(
                    SideAxisSpec(
                        column=col,
                        axis_type=spec.get("axis_type", "qualitative"),
                        palette=spec.get("palette"),
                    )
                )

        plot_config = PlotConfig(
            region=region,
            samples=sample_list,
            sort_by=sort_by_list or [],
            side_axes=sa_specs,
            colormap=colormap,
            panel_height_mode=panel_height_mode,
            max_reads=max_reads,
            min_mapq=min_mapq,
            min_cpg_coverage=min_cpg_coverage,
            min_cpgs_per_read=min_cpgs_per_read,
            nan_weight=nan_weight,
            gtf=gtf,
            gene_types=gene_types_list or ["protein_coding"],
            figsize=figsize_tuple,
            output=output,
            output_format=output_format,
            dpi=dpi,
        )
        plot_config.to_yaml(dump_config)
        click.echo(f"Configuration dumped to: {dump_config}")


if __name__ == "__main__":
    main()
