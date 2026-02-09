"""
Tests for utility functions.
"""

import pytest
from pathlib import Path
from cpgplotter.utils.io import parse_bam_spec, parse_figsize
from cpgplotter.utils.validation import validate_region


class TestIOUtils:
    """Test suite for I/O utility functions."""

    def test_parse_bam_spec_with_label(self):
        """Test parsing BAM spec with explicit label."""
        label, path = parse_bam_spec("Tumor:/data/tumor.bam")
        assert label == "Tumor"
        assert path == Path("/data/tumor.bam")

    def test_parse_bam_spec_without_label(self):
        """Test parsing BAM spec without label (derive from filename)."""
        label, path = parse_bam_spec("/data/sample.sorted.bam")
        assert label == "sample"
        assert path == Path("/data/sample.sorted.bam")

    def test_parse_figsize_valid(self):
        """Test parsing valid figsize string."""
        width, height = parse_figsize("14,8")
        assert width == 14.0
        assert height == 8.0

    def test_parse_figsize_invalid(self):
        """Test parsing invalid figsize string."""
        with pytest.raises(ValueError):
            parse_figsize("14")

        with pytest.raises(ValueError):
            parse_figsize("invalid")


class TestValidation:
    """Test suite for validation functions."""

    def test_validate_region_valid(self):
        """Test validation of valid region string."""
        chrom, start, end = validate_region("chr7:1072064-1101499")
        assert chrom == "chr7"
        assert start == 1072064
        assert end == 1101499

    def test_validate_region_invalid_format(self):
        """Test validation of invalid region format."""
        with pytest.raises(ValueError):
            validate_region("chr7_1072064-1101499")  # Wrong separator

        with pytest.raises(ValueError):
            validate_region("chr7:1072064")  # Missing end

    def test_validate_region_invalid_coords(self):
        """Test validation of invalid coordinates."""
        with pytest.raises(ValueError):
            validate_region("chr7:1101499-1072064")  # start > end
