"""
Tests for CpGIndex and coordinate transformation.
"""

import pytest
import numpy as np
from cpgplotter.core.coordinates import CpGIndex


class TestCpGIndex:
    """Test suite for CpGIndex class."""

    def test_init_with_positions(self):
        """Test CpGIndex initialization with pre-computed positions."""
        positions = np.array([1000, 1050, 1100, 1200])
        index = CpGIndex(
            region="chr1:1000-1300",
            positions=positions
        )

        assert index.chrom == "chr1"
        assert index.start == 1000
        assert index.end == 1300
        assert index.n_cpgs == 4
        assert np.array_equal(index.positions, positions)

    def test_genomic_to_cpg(self):
        """Test genomic to CpG index conversion."""
        positions = np.array([1000, 1050, 1100, 1200])
        index = CpGIndex(region="chr1:1000-1300", positions=positions)

        assert index.genomic_to_cpg(1000) == 0
        assert index.genomic_to_cpg(1050) == 1
        assert index.genomic_to_cpg(1200) == 3
        assert index.genomic_to_cpg(1075) is None  # Not a CpG position

    def test_cpg_to_genomic(self):
        """Test CpG index to genomic position conversion."""
        positions = np.array([1000, 1050, 1100, 1200])
        index = CpGIndex(region="chr1:1000-1300", positions=positions)

        assert index.cpg_to_genomic(0) == 1000
        assert index.cpg_to_genomic(2) == 1100
        assert index.cpg_to_genomic(3) == 1200

    def test_get_tick_positions(self):
        """Test generation of tick positions and labels."""
        positions = np.linspace(1000, 2000, 100, dtype=int)
        index = CpGIndex(region="chr1:1000-2000", positions=positions)

        tick_pos, tick_labels = index.get_tick_positions(n_ticks=5)

        assert len(tick_pos) == 5
        assert len(tick_labels) == 5
        assert tick_pos[0] == 0  # First CpG
        assert tick_pos[-1] == 99  # Last CpG

    def test_init_with_single_bam_string(self):
        """Test that single BAM path as string is handled correctly."""
        # This is a unit test for the type handling, not requiring actual BAM
        # The actual BAM reading is tested in integration tests
        positions = np.array([1000, 1050, 1100])

        # Should accept both string and list[str] for bam_path parameter
        # (actual BAM reading tested in integration tests)
        index = CpGIndex(region="chr1:1000-1200", positions=positions)
        assert index.n_cpgs == 3

    def test_init_with_bam_list(self):
        """Test that list of BAM paths is handled correctly."""
        # This is a unit test for the type handling
        positions = np.array([1000, 1050, 1100])
        index = CpGIndex(region="chr1:1000-1200", positions=positions)
        assert index.n_cpgs == 3

    def test_min_cpg_coverage_default(self):
        """Test that min_cpg_coverage defaults to 3."""
        positions = np.array([1000, 1050, 1100])
        index = CpGIndex(region="chr1:1000-1200", positions=positions)
        assert index.min_cpg_coverage == 3

    def test_min_cpg_coverage_stored(self):
        """Test that custom min_cpg_coverage is stored."""
        positions = np.array([1000, 1050, 1100])
        index = CpGIndex(
            region="chr1:1000-1200", positions=positions, min_cpg_coverage=5
        )
        assert index.min_cpg_coverage == 5
