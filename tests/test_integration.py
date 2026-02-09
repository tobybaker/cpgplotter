"""
Integration tests using real BAM data.

These tests use the test_data/COLO829.GRCh38_filtered.bam file to test
the complete pipeline with real methylation data.
"""

import pytest
import random
from pathlib import Path
import numpy as np
import polars as pl
from cpgplotter import plot_methylation, CpGIndex, extract_methylation
from cpgplotter.core.annotations import ReadAnnotations, ReadIntervals

# Test data paths
TEST_DATA_DIR = Path(__file__).parent.parent / "test_data"
TEST_BAM = TEST_DATA_DIR / "COLO829.GRCh38_filtered.bam"
TEST_REGIONS = TEST_DATA_DIR / "regions.tsv"
TEST_READ_ANNOTATIONS = TEST_DATA_DIR / "read_annotations.tsv"
TEST_READ_REGIONS = TEST_DATA_DIR / "read_regions.tsv"


def get_random_40kb_regions(n_regions: int = 2, seed: int = 42) -> list[str]:
    """
    Generate random 40kb sub-regions from the test regions file.

    Args:
        n_regions: Number of random 40kb regions to generate
        seed: Random seed for reproducibility

    Returns:
        List of region strings in format "chr:start-end"
    """
    random.seed(seed)

    # Read available regions
    regions_df = pl.read_csv(TEST_REGIONS, separator="\t")

    # Sample a few source regions
    sampled_regions = regions_df.sample(n=min(n_regions, len(regions_df)), seed=seed)

    # Generate 40kb sub-regions within each sampled region
    result_regions = []
    for row in sampled_regions.iter_rows(named=True):
        chrom = row["chrom"]
        region_start = row["start"]
        region_end = row["end"]

        # Generate random 40kb window within this region
        # Ensure we have at least 40kb of space
        if region_end - region_start < 40000:
            continue

        # Random start position (leaving room for 40kb window)
        max_start = region_end - 40000
        start = random.randint(region_start, max_start)
        end = start + 40000

        result_regions.append(f"{chrom}:{start}-{end}")

    return result_regions


class TestIntegration:
    """Integration tests with real BAM data."""

    @pytest.fixture(scope="class")
    def test_regions(self):
        """Generate test regions for the test suite."""
        return get_random_40kb_regions(n_regions=3, seed=42)

    def test_bam_file_exists(self):
        """Verify test BAM file exists and is accessible."""
        assert TEST_BAM.exists(), f"Test BAM not found at {TEST_BAM}"
        assert TEST_BAM.with_suffix(".bam.bai").exists(), "BAM index not found"

    def test_regions_file_exists(self):
        """Verify regions file exists and is readable."""
        assert TEST_REGIONS.exists(), f"Regions file not found at {TEST_REGIONS}"
        df = pl.read_csv(TEST_REGIONS, separator="\t")
        assert len(df) > 0, "Regions file is empty"
        assert all(col in df.columns for col in ["chrom", "start", "end"])

    def test_annotation_files_exist(self):
        """Verify annotation files exist and are readable."""
        assert TEST_READ_ANNOTATIONS.exists(), f"Read annotations not found at {TEST_READ_ANNOTATIONS}"
        assert TEST_READ_REGIONS.exists(), f"Read regions not found at {TEST_READ_REGIONS}"

        # Verify read annotations format
        annot_df = pl.read_csv(TEST_READ_ANNOTATIONS, separator="\t")
        assert len(annot_df) > 0, "Read annotations file is empty"
        assert "read_name" in annot_df.columns

        # Verify read regions format
        regions_df = pl.read_csv(TEST_READ_REGIONS, separator="\t")
        assert len(regions_df) > 0, "Read regions file is empty"
        assert all(col in regions_df.columns for col in ["read_name", "start", "end"])

    def test_random_region_generation(self):
        """Test that random region generation works correctly."""
        regions = get_random_40kb_regions(n_regions=5, seed=123)

        assert len(regions) > 0, "No regions generated"

        # Verify each region is properly formatted and is 40kb
        for region in regions:
            assert ":" in region, f"Invalid region format: {region}"
            assert "-" in region, f"Invalid region format: {region}"

            chrom, coords = region.split(":")
            start, end = map(int, coords.split("-"))

            assert end - start == 40000, f"Region {region} is not 40kb"
            assert start >= 0, f"Invalid start position in {region}"

    def test_cpg_index_from_bam(self, test_regions):
        """Test CpGIndex construction from BAM file."""
        region = test_regions[0]

        cpg_index = CpGIndex(region=region, bam_path=str(TEST_BAM))

        assert cpg_index.n_cpgs > 0, "No CpGs found in region"
        assert cpg_index.chrom in region, "Chromosome mismatch"

    def test_methylation_extraction(self, test_regions):
        """Test methylation data extraction from BAM."""
        region = test_regions[0]

        # First create CpG index
        cpg_index = CpGIndex(region=region, bam_path=str(TEST_BAM))

        # Extract methylation data
        meth_data = extract_methylation(
            bam_path=str(TEST_BAM),
            cpg_index=cpg_index,
            min_mapq=0
        )

        assert meth_data.n_reads > 0, "No reads found in region"
        assert meth_data.n_cpgs == cpg_index.n_cpgs, "CpG count mismatch"
        assert meth_data.matrix.shape == (meth_data.n_reads, meth_data.n_cpgs)

    def test_cpg_index_with_multiple_bams(self, test_regions):
        """Test CpGIndex construction from multiple BAM files (shared coordinate system)."""
        region = test_regions[0]

        # Create index from single BAM
        cpg_index_single = CpGIndex(region=region, bam_path=str(TEST_BAM))

        # Create index from list with single BAM (should be equivalent)
        cpg_index_list = CpGIndex(region=region, bam_path=[str(TEST_BAM)])

        # Both should produce the same CpG positions
        assert cpg_index_single.n_cpgs == cpg_index_list.n_cpgs
        assert np.array_equal(cpg_index_single.positions, cpg_index_list.positions)

        # In a real multi-sample scenario, you might have:
        # cpg_index_shared = CpGIndex(region=region, bam_path=[tumor_bam, normal_bam])
        # This would create a shared coordinate system with CpGs from both samples

    def test_min_cpg_coverage_filters_spurious_sites(self, test_regions):
        """Test that min_cpg_coverage filters out low-coverage CpG positions."""
        region = test_regions[0]

        # No filtering (min_cpg_coverage=1 keeps everything)
        index_all = CpGIndex(region=region, bam_path=str(TEST_BAM), min_cpg_coverage=1)

        # Default filtering (min_cpg_coverage=3)
        index_filtered = CpGIndex(region=region, bam_path=str(TEST_BAM), min_cpg_coverage=3)

        # Filtered index should have fewer or equal CpGs
        assert index_filtered.n_cpgs <= index_all.n_cpgs
        # All filtered positions should be a subset of unfiltered positions
        assert set(index_filtered.positions).issubset(set(index_all.positions))

    def test_min_cpgs_per_read_default(self, test_regions):
        """Test that reads with fewer than 5 CpGs (default) are excluded."""
        region = test_regions[0]

        cpg_index = CpGIndex(region=region, bam_path=str(TEST_BAM))
        meth_data = extract_methylation(bam_path=str(TEST_BAM), cpg_index=cpg_index)

        # Every row should have at least min_cpgs_per_read (default=5) non-NaN values
        cpgs_per_read = np.sum(~np.isnan(meth_data.matrix), axis=1)
        assert np.all(cpgs_per_read >= 5), "Found reads with fewer than 5 CpGs"

    def test_min_cpgs_per_read_configurable(self, test_regions):
        """Test that min_cpgs_per_read threshold is respected at different values."""
        region = test_regions[0]
        cpg_index = CpGIndex(region=region, bam_path=str(TEST_BAM))

        # Permissive: keep reads with >= 1 CpG
        meth_loose = extract_methylation(
            bam_path=str(TEST_BAM), cpg_index=cpg_index, min_cpgs_per_read=1,
        )
        # Strict: keep reads with >= 10 CpGs
        meth_strict = extract_methylation(
            bam_path=str(TEST_BAM), cpg_index=cpg_index, min_cpgs_per_read=10,
        )

        # Stricter threshold should keep fewer or equal reads
        assert meth_strict.n_reads <= meth_loose.n_reads

        # Verify each matrix respects its threshold
        if meth_loose.n_reads > 0:
            assert np.all(np.sum(~np.isnan(meth_loose.matrix), axis=1) >= 1)
        if meth_strict.n_reads > 0:
            assert np.all(np.sum(~np.isnan(meth_strict.matrix), axis=1) >= 10)

    @pytest.mark.skip(reason="Full pipeline not yet implemented")
    def test_plot_single_region(self, test_regions, tmp_path):
        """Test plotting a single 40kb region."""
        region = test_regions[0]
        output_file = tmp_path / "test_plot.png"

        fig = plot_methylation(
            region=region,
            bams=str(TEST_BAM),
            output=str(output_file),
            min_mapq=10
        )

        assert output_file.exists(), "Output file not created"
        assert fig is not None, "Figure not returned"

    @pytest.mark.skip(reason="Full pipeline not yet implemented")
    def test_plot_multiple_regions(self, test_regions, tmp_path):
        """Test plotting multiple 40kb regions from different chromosomes."""
        for i, region in enumerate(test_regions[:2]):
            output_file = tmp_path / f"test_plot_region_{i}.png"

            fig = plot_methylation(
                region=region,
                bams={f"COLO829_region{i}": str(TEST_BAM)},
                output=str(output_file),
                min_mapq=10,
                colormap="RdBu_r"
            )

            assert output_file.exists(), f"Output file {i} not created"

    def test_load_annotations_for_region(self, test_regions):
        """Test loading annotations for reads in a specific region."""
        region = test_regions[0]

        # Load all annotations
        all_annotations = ReadAnnotations.from_tsv(TEST_READ_ANNOTATIONS)

        # Get reads from BAM in this region
        import pysam
        bam = pysam.AlignmentFile(str(TEST_BAM), "rb")
        region_reads = {read.query_name for read in bam.fetch(region=region) if not read.is_unmapped}
        bam.close()

        # Filter annotations to this region
        region_annotations = all_annotations.filter_reads(region_reads)

        assert len(region_annotations.data) > 0, "No annotations found for region reads"
        assert len(region_annotations.data) <= len(all_annotations.data)
        assert all(name in region_reads for name in region_annotations.data["read_name"])

    def test_load_intervals_for_region(self, test_regions):
        """Test loading intervals for reads in a specific region."""
        region = test_regions[0]
        chrom, coords = region.split(":")
        start, end = map(int, coords.split("-"))

        # Load all intervals
        all_intervals = ReadIntervals.from_tsv(TEST_READ_REGIONS)

        # Get reads from BAM in this region
        import pysam
        bam = pysam.AlignmentFile(str(TEST_BAM), "rb")
        region_reads = {read.query_name for read in bam.fetch(region=region) if not read.is_unmapped}
        bam.close()

        # Filter intervals to this region
        region_intervals = all_intervals.filter_reads(region_reads)

        assert len(region_intervals.data) > 0, "No intervals found for region reads"

        # Verify intervals exist for the read names (they may not overlap the specific
        # query window since intervals span the full read, which may extend beyond the window)
        interval_read_names = set(region_intervals.data["read_name"].unique())
        assert interval_read_names.issubset(region_reads), "Intervals contain unexpected reads"

        # Verify each read has the expected number of intervals (10)
        intervals_per_read = region_intervals.data.group_by("read_name").len()
        assert all(intervals_per_read["len"] == 10), "Not all reads have exactly 10 intervals"

    def test_annotations_and_intervals_consistency(self, test_regions):
        """Test that annotations and intervals are consistent for the same reads."""
        region = test_regions[0]

        # Load both annotations and intervals
        annotations = ReadAnnotations.from_tsv(TEST_READ_ANNOTATIONS)
        intervals = ReadIntervals.from_tsv(TEST_READ_REGIONS)

        # Get reads from BAM
        import pysam
        bam = pysam.AlignmentFile(str(TEST_BAM), "rb")
        region_reads = {read.query_name for read in bam.fetch(region=region) if not read.is_unmapped}
        bam.close()

        # Filter both to region reads
        region_annotations = annotations.filter_reads(region_reads)
        region_intervals = intervals.filter_reads(region_reads)

        # Get unique reads from each
        annot_reads = set(region_annotations.data["read_name"])
        interval_reads = set(region_intervals.data["read_name"].unique())

        # All reads with intervals should have annotations
        assert interval_reads.issubset(annot_reads), "Some reads have intervals but no annotations"

        # Verify we can cross-reference
        common_reads = annot_reads & interval_reads
        assert len(common_reads) > 0, "No reads have both annotations and intervals"


def test_print_test_regions(capsys):
    """Print out the test regions for manual inspection."""
    regions = get_random_40kb_regions(n_regions=5, seed=42)

    print("\nGenerated test regions (40kb each):")
    for i, region in enumerate(regions, 1):
        print(f"  {i}. {region}")

    captured = capsys.readouterr()
    assert "Generated test regions" in captured.out


if __name__ == "__main__":
    # When run directly, print some test regions
    print("Test data location:", TEST_DATA_DIR)
    print("BAM file:", TEST_BAM)
    print("BAM exists:", TEST_BAM.exists())
    print("\nRandom 40kb test regions:")

    regions = get_random_40kb_regions(n_regions=5, seed=42)
    for i, region in enumerate(regions, 1):
        print(f"  {i}. {region}")
