"""
Per-read annotation management.

This module implements ReadAnnotations and ReadIntervals classes for handling
both scalar (per-read) and interval (sub-read region) annotations.
"""

from pathlib import Path
from typing import Optional, Literal
import polars as pl
import pysam


class ReadAnnotations:
    """
    Container for per-read scalar annotations.

    ReadAnnotations stores qualitative and quantitative annotations for reads,
    such as haplotype assignment, mapping quality, or cluster membership.
    Annotations can be loaded from TSV files or extracted from BAM tags.

    Attributes:
        data: DataFrame indexed by read_name with annotation columns
        column_types: Mapping of column names to annotation types
    """

    def __init__(
        self,
        data: pl.DataFrame,
        column_types: Optional[dict[str, Literal["qualitative", "quantitative"]]] = None,
    ):
        """
        Initialize ReadAnnotations from a DataFrame.

        Args:
            data: DataFrame with 'read_name' column and annotation columns
            column_types: Optional explicit type mapping for columns
                         If None, types are inferred from dtypes
        """
        if "read_name" not in data.columns:
            raise ValueError("DataFrame must have 'read_name' column")

        self.data = data

        # Infer column types if not provided
        if column_types is None:
            self.column_types = self._infer_column_types()
        else:
            self.column_types = column_types

    def _infer_column_types(self) -> dict[str, Literal["qualitative", "quantitative"]]:
        """
        Infer annotation types from DataFrame dtypes.

        Returns:
            Mapping of column names to inferred types
        """
        types = {}
        for col in self.data.columns:
            if col == "read_name":
                continue

            dtype = self.data[col].dtype
            if dtype in [pl.Float32, pl.Float64, pl.Int32, pl.Int64, pl.UInt32, pl.UInt64]:
                types[col] = "quantitative"
            else:
                types[col] = "qualitative"

        return types

    @classmethod
    def from_tsv(cls, path: Path) -> "ReadAnnotations":
        """
        Load annotations from TSV file.

        Expected format:
            read_name    annotation1    annotation2    ...
            read001      value1         value2         ...

        Args:
            path: Path to TSV file

        Returns:
            ReadAnnotations instance
        """
        data = pl.read_csv(path, separator="\t")
        return cls(data)

    @classmethod
    def from_bam_tags(
        cls,
        bam_path: str,
        region: str,
        tags: dict[str, str],
    ) -> "ReadAnnotations":
        """
        Extract annotations from BAM tags.

        Args:
            bam_path: Path to BAM file
            region: Genomic region to extract reads from
            tags: Mapping of annotation_name -> tag_name
                 Special value "mapq" extracts mapping quality

        Returns:
            ReadAnnotations instance

        Example:
            >>> annot = ReadAnnotations.from_bam_tags(
            ...     "sample.bam",
            ...     "chr7:1072064-1101499",
            ...     {"haplotype": "HP", "mapq": "mapq"}
            ... )
        """
        chrom, coords = region.split(":")
        start, end = map(int, coords.split("-"))

        rows = []
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(chrom, start, end):
                if read.is_unmapped:
                    continue

                row = {"read_name": read.query_name}
                for annot_name, source in tags.items():
                    if source.lower() == "mapq":
                        row[annot_name] = read.mapping_quality
                    else:
                        try:
                            row[annot_name] = read.get_tag(source)
                        except KeyError:
                            row[annot_name] = None
                rows.append(row)

        if not rows:
            schema = {"read_name": pl.Utf8}
            for annot_name in tags:
                schema[annot_name] = pl.Utf8
            return cls(pl.DataFrame(schema=schema))

        data = pl.DataFrame(rows)
        return cls(data)

    def filter_reads(self, read_names: set[str]) -> "ReadAnnotations":
        """
        Filter annotations to specified read names.

        Args:
            read_names: Set of read names to keep

        Returns:
            New ReadAnnotations with filtered data
        """
        filtered_data = self.data.filter(pl.col("read_name").is_in(read_names))
        return ReadAnnotations(filtered_data, self.column_types)

    def get_column(self, column: str) -> pl.Series:
        """Get a specific annotation column."""
        return self.data[column]

    def __repr__(self) -> str:
        n_reads = len(self.data)
        n_cols = len(self.data.columns) - 1  # Exclude read_name
        return f"ReadAnnotations(n_reads={n_reads}, n_columns={n_cols})"


class ReadIntervals:
    """
    Container for per-read interval annotations.

    ReadIntervals stores sub-read regions of interest (e.g., FIRE peaks,
    accessible chromatin regions) in both genomic and CpG coordinate spaces.

    Attributes:
        data: DataFrame with columns: read_name, start, end, label, score
              start/end are in genomic coordinates
        cpg_data: Optional DataFrame with CpG space coordinates after transformation
    """

    def __init__(self, data: pl.DataFrame):
        """
        Initialize ReadIntervals from a DataFrame.

        Args:
            data: DataFrame with required columns: read_name, start, end
                 Optional columns: label, score
        """
        required_cols = {"read_name", "start", "end"}
        if not required_cols.issubset(set(data.columns)):
            raise ValueError(f"DataFrame must have columns: {required_cols}")

        self.data = data
        self.cpg_data: Optional[pl.DataFrame] = None

    @classmethod
    def from_tsv(cls, path: Path) -> "ReadIntervals":
        """
        Load intervals from TSV file with header.

        Expected format (tab-separated with header):
            read_name    start    end    [label]    [score]

        Args:
            path: Path to TSV file

        Returns:
            ReadIntervals instance
        """
        data = pl.read_csv(path, separator="\t")
        return cls(data)

    @classmethod
    def from_bed(cls, path: Path) -> "ReadIntervals":
        """
        Load intervals from BED file (no header).

        Expected format (tab-separated, no header):
            read_name    start    end    [label]    [score]

        Args:
            path: Path to BED file

        Returns:
            ReadIntervals instance
        """
        # Read BED file - handle both 3-column and 5-column formats
        data = pl.read_csv(path, separator="\t", has_header=False)

        # Rename columns based on number of columns present
        n_cols = len(data.columns)
        col_names = ["read_name", "start", "end", "label", "score"][:n_cols]

        # Create rename mapping from auto-generated column names
        rename_map = {f"column_{i+1}": name for i, name in enumerate(col_names)}
        data = data.rename(rename_map)

        return cls(data)

    def transform_to_cpg_space(self, cpg_index) -> None:
        """
        Transform intervals from genomic to CpG coordinate space.

        This modifies the object in-place by populating the cpg_data attribute.

        Args:
            cpg_index: CpGIndex object for coordinate transformation
        """
        self.cpg_data = cpg_index.transform_intervals(self.data)

    def filter_reads(self, read_names: set[str]) -> "ReadIntervals":
        """
        Filter intervals to specified read names.

        Args:
            read_names: Set of read names to keep

        Returns:
            New ReadIntervals with filtered data
        """
        filtered_data = self.data.filter(pl.col("read_name").is_in(read_names))
        result = ReadIntervals(filtered_data)
        if self.cpg_data is not None:
            result.cpg_data = self.cpg_data.filter(pl.col("read_name").is_in(read_names))
        return result

    def __repr__(self) -> str:
        n_intervals = len(self.data)
        n_reads = self.data["read_name"].n_unique()
        return f"ReadIntervals(n_intervals={n_intervals}, n_reads={n_reads})"
