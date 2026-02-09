"""
Tests for annotation handling with real test data.

These tests use the generated test annotation files to validate
ReadAnnotations and ReadIntervals functionality.
"""

import pytest
from pathlib import Path
import polars as pl
import numpy as np
from cpgplotter.core.annotations import ReadAnnotations, ReadIntervals

# Test data paths
TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
TEST_BAM = TEST_DATA_DIR / "COLO829.GRCh38_filtered.bam"
TEST_READ_ANNOTATIONS = TEST_DATA_DIR / "read_annotations.tsv"
TEST_READ_REGIONS = TEST_DATA_DIR / "read_regions.tsv"


class TestReadAnnotations:
    """Test suite for ReadAnnotations with real data."""

    @pytest.fixture(scope="class")
    def annotations_df(self):
        """Load test annotations DataFrame."""
        return pl.read_csv(TEST_READ_ANNOTATIONS, separator="\t")

    @pytest.fixture(scope="class")
    def annotations(self, annotations_df):
        """Create ReadAnnotations from test data."""
        return ReadAnnotations(annotations_df)

    def test_annotation_files_exist(self):
        """Verify test annotation files exist."""
        assert TEST_READ_ANNOTATIONS.exists(), f"Annotations not found at {TEST_READ_ANNOTATIONS}"
        assert TEST_READ_REGIONS.exists(), f"Regions not found at {TEST_READ_REGIONS}"

    def test_load_from_tsv(self):
        """Test loading annotations from TSV file."""
        annotations = ReadAnnotations.from_tsv(TEST_READ_ANNOTATIONS)

        assert annotations.data is not None
        assert "read_name" in annotations.data.columns
        assert len(annotations.data) > 0

    def test_annotation_schema(self, annotations_df):
        """Test that annotation file has expected schema."""
        expected_columns = {"read_name", "binary_label", "qualitative_label", "quantitative_label"}
        assert expected_columns.issubset(set(annotations_df.columns))

    def test_binary_label_values(self, annotations_df):
        """Test binary label contains only 0 and 1."""
        unique_values = annotations_df["binary_label"].unique().sort()
        assert len(unique_values) == 2
        assert all(v in [0, 1] for v in unique_values)

    def test_qualitative_label_categories(self, annotations_df):
        """Test qualitative label has expected categories."""
        unique_categories = set(annotations_df["qualitative_label"].unique())
        expected_categories = {"A", "B", "C", "D", "E"}
        assert unique_categories == expected_categories

    def test_quantitative_label_range(self, annotations_df):
        """Test quantitative label is in [0, 1] range."""
        values = annotations_df["quantitative_label"].to_numpy()
        assert values.min() >= 0.0
        assert values.max() <= 1.0

    def test_quantitative_label_distribution(self, annotations_df):
        """Test quantitative label follows uniform distribution."""
        values = annotations_df["quantitative_label"].to_numpy()

        # For uniform distribution, mean should be around 0.5
        mean = values.mean()
        assert 0.45 < mean < 0.55, f"Mean {mean} not close to 0.5"

        # Standard deviation should be around 0.289
        std = values.std()
        assert 0.25 < std < 0.32, f"Std {std} not close to 0.289"

    def test_column_type_inference(self, annotations):
        """Test automatic type inference for annotation columns."""
        # Binary label should be quantitative (it's numeric)
        assert annotations.column_types["binary_label"] == "quantitative"

        # Qualitative label should be qualitative (it's string)
        assert annotations.column_types["qualitative_label"] == "qualitative"

        # Quantitative label should be quantitative (it's float)
        assert annotations.column_types["quantitative_label"] == "quantitative"

    def test_get_column(self, annotations):
        """Test retrieving specific annotation columns."""
        haplotype_col = annotations.get_column("qualitative_label")
        assert isinstance(haplotype_col, pl.Series)
        assert len(haplotype_col) == len(annotations.data)

    def test_filter_reads(self, annotations):
        """Test filtering annotations to specific reads."""
        # Get first 100 read names
        sample_reads = set(annotations.data["read_name"].head(100))

        filtered = annotations.filter_reads(sample_reads)

        assert len(filtered.data) == len(sample_reads)
        assert all(name in sample_reads for name in filtered.data["read_name"])

    def test_filter_reads_preserves_types(self, annotations):
        """Test that filtering preserves column types."""
        sample_reads = set(annotations.data["read_name"].head(50))
        filtered = annotations.filter_reads(sample_reads)

        assert filtered.column_types == annotations.column_types

    def test_repr(self, annotations):
        """Test string representation."""
        repr_str = repr(annotations)
        assert "ReadAnnotations" in repr_str
        assert "n_reads=" in repr_str
        assert "n_columns=" in repr_str


class TestReadIntervals:
    """Test suite for ReadIntervals with real data."""

    @pytest.fixture(scope="class")
    def intervals_df(self):
        """Load test intervals DataFrame."""
        return pl.read_csv(TEST_READ_REGIONS, separator="\t")

    @pytest.fixture(scope="class")
    def intervals(self, intervals_df):
        """Create ReadIntervals from test data."""
        return ReadIntervals(intervals_df)

    def test_load_from_tsv(self):
        """Test loading intervals from TSV file."""
        intervals = ReadIntervals.from_tsv(TEST_READ_REGIONS)

        assert intervals.data is not None
        assert "read_name" in intervals.data.columns
        assert "start" in intervals.data.columns
        assert "end" in intervals.data.columns
        assert len(intervals.data) > 0

    def test_interval_schema(self, intervals_df):
        """Test that interval file has expected schema."""
        required_columns = {"read_name", "start", "end"}
        assert required_columns.issubset(set(intervals_df.columns))

        # Optional columns
        assert "label" in intervals_df.columns
        assert "score" in intervals_df.columns

    def test_interval_coordinates_valid(self, intervals_df):
        """Test that interval coordinates are valid."""
        # Start should be less than end for all intervals
        assert all(intervals_df["start"] < intervals_df["end"])

        # Coordinates should be positive
        assert intervals_df["start"].min() >= 0
        assert intervals_df["end"].min() > 0

    def test_interval_lengths(self, intervals_df):
        """Test interval length distribution."""
        lengths = intervals_df["end"] - intervals_df["start"]

        # Most intervals should be between 100-1000bp (as generated)
        # Allow some shorter ones due to read boundary constraints
        assert lengths.min() >= 0
        assert lengths.max() <= 1000

        # Mean should be around 500bp
        mean_length = lengths.mean()
        assert 400 < mean_length < 650, f"Mean length {mean_length} outside expected range"

    def test_intervals_per_read(self, intervals_df):
        """Test that reads have expected number of intervals."""
        intervals_per_read = intervals_df.group_by("read_name").len()

        # Each read should have exactly 10 regions (5 with label A, 5 with label B)
        counts = intervals_per_read["len"]
        assert counts.min() == 10
        assert counts.max() == 10

        # Mean should be exactly 10
        mean_count = counts.mean()
        assert mean_count == 10.0, f"Mean count {mean_count} is not 10.0"

    def test_interval_labels(self, intervals_df):
        """Test interval label categories."""
        unique_labels = set(intervals_df["label"].unique())
        expected_labels = {"A", "B"}
        assert unique_labels == expected_labels

        # Each read should have 5 of each label
        label_counts = intervals_df.group_by(["read_name", "label"]).len()
        # All counts should be 5
        assert all(label_counts["len"] == 5)

    def test_interval_scores(self, intervals_df):
        """Test interval scores are in valid range."""
        scores = intervals_df["score"].to_numpy()

        # Scores should be between 0.5 and 1.0 (as generated)
        assert scores.min() >= 0.5
        assert scores.max() <= 1.0

    def test_filter_reads(self, intervals):
        """Test filtering intervals to specific reads."""
        # Get first 50 unique read names
        unique_reads = intervals.data["read_name"].unique()
        sample_reads = set(unique_reads.head(50))

        filtered = intervals.filter_reads(sample_reads)

        assert len(filtered.data) > 0
        assert all(name in sample_reads for name in filtered.data["read_name"])

    def test_filter_reads_preserves_intervals(self, intervals):
        """Test that filtering preserves interval structure."""
        unique_reads = intervals.data["read_name"].unique()
        sample_reads = set(unique_reads.head(10))

        original_counts = (
            intervals.data.filter(pl.col("read_name").is_in(sample_reads))
            .group_by("read_name")
            .len()
        )

        filtered = intervals.filter_reads(sample_reads)
        filtered_counts = filtered.data.group_by("read_name").len()

        # Should have same number of intervals per read
        assert len(original_counts) == len(filtered_counts)

    def test_repr(self, intervals):
        """Test string representation."""
        repr_str = repr(intervals)
        assert "ReadIntervals" in repr_str
        assert "n_intervals=" in repr_str
        assert "n_reads=" in repr_str


class TestAnnotationIntegration:
    """Integration tests combining annotations and intervals."""

    @pytest.fixture(scope="class")
    def annotations(self):
        """Load read annotations."""
        return ReadAnnotations.from_tsv(TEST_READ_ANNOTATIONS)

    @pytest.fixture(scope="class")
    def intervals(self):
        """Load read intervals."""
        return ReadIntervals.from_tsv(TEST_READ_REGIONS)

    def test_matching_read_names(self, annotations, intervals):
        """Test that annotations and intervals have matching reads."""
        annot_reads = set(annotations.data["read_name"])
        interval_reads = set(intervals.data["read_name"].unique())

        # All reads with intervals should have annotations
        assert interval_reads.issubset(annot_reads)

    def test_combined_filtering(self, annotations, intervals):
        """Test filtering both annotations and intervals to same read set."""
        # Get reads that appear in both
        annot_reads = set(annotations.data["read_name"])
        interval_reads = set(intervals.data["read_name"].unique())
        common_reads = annot_reads & interval_reads

        # Sample 100 reads
        sample_reads = set(list(common_reads)[:100])

        filtered_annot = annotations.filter_reads(sample_reads)
        filtered_intervals = intervals.filter_reads(sample_reads)

        # Both should have the same read names
        filtered_annot_reads = set(filtered_annot.data["read_name"])
        filtered_interval_reads = set(filtered_intervals.data["read_name"].unique())

        assert filtered_annot_reads == sample_reads
        assert filtered_interval_reads.issubset(sample_reads)

    def test_stratify_by_annotation(self, annotations, intervals):
        """Test stratifying intervals by annotation value."""
        # Get intervals for reads with qualitative_label == "A"
        reads_label_a = set(
            annotations.data.filter(pl.col("qualitative_label") == "A")["read_name"]
        )

        intervals_label_a = intervals.filter_reads(reads_label_a)

        assert len(intervals_label_a.data) > 0
        assert all(
            name in reads_label_a for name in intervals_label_a.data["read_name"]
        )

    def test_quantitative_annotation_correlation_with_intervals(
        self, annotations, intervals
    ):
        """Test relationship between quantitative annotations and interval counts."""
        # Count intervals per read
        interval_counts = intervals.data.group_by("read_name").len()

        # Join with quantitative annotation
        merged = interval_counts.join(
            annotations.data.select(["read_name", "quantitative_label"]),
            on="read_name",
            how="inner",
        )

        # Should have data for analysis
        assert len(merged) > 0
        assert "len" in merged.columns
        assert "quantitative_label" in merged.columns

        # The relationship should be random (no correlation expected)
        # But we can verify both columns have reasonable ranges
        assert merged["len"].min() == 10
        assert merged["len"].max() == 10
        assert merged["quantitative_label"].min() >= 0.0
        assert merged["quantitative_label"].max() <= 1.0


class TestAnnotationEdgeCases:
    """Test edge cases and error handling."""

    def test_missing_read_name_column(self):
        """Test error when DataFrame lacks read_name column."""
        bad_df = pl.DataFrame({"foo": [1, 2, 3], "bar": [4, 5, 6]})

        with pytest.raises(ValueError, match="read_name"):
            ReadAnnotations(bad_df)

    def test_missing_required_interval_columns(self):
        """Test error when DataFrame lacks required interval columns."""
        bad_df = pl.DataFrame({"read_name": ["r1", "r2"], "start": [100, 200]})

        with pytest.raises(ValueError, match="DataFrame must have columns"):
            ReadIntervals(bad_df)

    def test_empty_annotations(self):
        """Test handling of empty annotation DataFrame."""
        empty_df = pl.DataFrame(
            {"read_name": [], "annotation": []}, schema={"read_name": pl.Utf8, "annotation": pl.Int64}
        )

        annotations = ReadAnnotations(empty_df)
        assert len(annotations.data) == 0

    def test_empty_intervals(self):
        """Test handling of empty intervals DataFrame."""
        empty_df = pl.DataFrame(
            {"read_name": [], "start": [], "end": []},
            schema={"read_name": pl.Utf8, "start": pl.Int64, "end": pl.Int64},
        )

        intervals = ReadIntervals(empty_df)
        assert len(intervals.data) == 0

    def test_filter_with_no_matches(self):
        """Test filtering with read names that don't exist."""
        annotations = ReadAnnotations.from_tsv(TEST_READ_ANNOTATIONS)
        fake_reads = {"fake_read_1", "fake_read_2", "fake_read_3"}

        filtered = annotations.filter_reads(fake_reads)
        assert len(filtered.data) == 0
