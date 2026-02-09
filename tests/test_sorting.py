"""
Tests for read sorting and clustering algorithms.
"""

import numpy as np
import polars as pl
import pytest
from scipy.cluster.hierarchy import linkage

from cpgplotter.core.annotations import ReadAnnotations
from cpgplotter.processing.sorting import (
    compute_hellinger_distance_matrix,
    hierarchical_cluster_reads,
    sort_reads,
)


class TestComputeHellingerDistanceMatrix:
    """Tests for pairwise Hellinger distance matrix computation."""

    def test_distance_matrix_small_overlap(self):
        """Test with small matrix with known overlaps."""
        matrix = np.array(
            [
                [0.0, 0.5, np.nan, np.nan],  # read 0
                [0.1, 0.6, np.nan, np.nan],  # read 1 (overlap with 0 on CpGs 0,1)
                [np.nan, np.nan, 0.8, 0.9],  # read 2 (no overlap with 0,1)
            ]
        )

        distances = compute_hellinger_distance_matrix(matrix)

        # Should have 3 distances: (0,1), (0,2), (1,2)
        assert len(distances) == 3

        # (0,1) should have small distance (overlap on CpGs 0,1)
        assert distances[0] < 1.0

        # (0,2) should be sqrt(2) (no overlap)
        assert distances[1] == pytest.approx(np.sqrt(2.0))

        # (1,2) should be sqrt(2) (no overlap)
        assert distances[2] == pytest.approx(np.sqrt(2.0))

    def test_distance_matrix_identical_reads(self):
        """Identical reads should have distance ~0."""
        matrix = np.array([[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]])

        distances = compute_hellinger_distance_matrix(matrix)

        assert len(distances) == 1
        assert distances[0] < 1e-10

    def test_distance_matrix_no_overlap(self):
        """All reads non-overlapping should give all distances = sqrt(2)."""
        matrix = np.array(
            [
                [0.5, np.nan, np.nan],
                [np.nan, 0.5, np.nan],
                [np.nan, np.nan, 0.5],
            ]
        )

        distances = compute_hellinger_distance_matrix(matrix)

        # 3 reads → 3 distances
        assert len(distances) == 3
        assert distances == pytest.approx(np.full(3, np.sqrt(2.0)))

    def test_distance_matrix_compatible_with_linkage(self):
        """Distance matrix should be compatible with scipy linkage."""
        # Create random matrix
        np.random.seed(42)
        matrix = np.random.rand(10, 20)

        distances = compute_hellinger_distance_matrix(matrix)

        # Should be able to create linkage without errors
        linkage_matrix = linkage(distances, method="ward")

        # Linkage matrix should have n-1 rows for n reads
        assert linkage_matrix.shape == (9, 4)

    def test_distance_matrix_single_read(self):
        """Single read should return empty distance array."""
        matrix = np.array([[0.5, 0.5, 0.5]])

        distances = compute_hellinger_distance_matrix(matrix)

        assert len(distances) == 0

    def test_distance_matrix_two_reads(self):
        """Two reads should return single distance."""
        matrix = np.array([[0.0, 0.5], [0.1, 0.6]])

        distances = compute_hellinger_distance_matrix(matrix)

        assert len(distances) == 1
        assert distances[0] > 0


class TestHierarchicalClusterReads:
    """Tests for hierarchical clustering integration."""

    def test_cluster_reads_basic(self):
        """Test basic clustering with overlap."""
        np.random.seed(42)
        matrix = np.random.rand(20, 50)

        indices = hierarchical_cluster_reads(matrix)

        # Should return all indices in some order
        assert len(indices) == 20
        assert set(indices) == set(range(20))

    def test_cluster_reads_subset(self):
        """Test clustering with read subset."""
        np.random.seed(42)
        matrix = np.random.rand(50, 100)

        # Cluster only first 10 reads
        subset_indices = np.arange(10)
        clustered = hierarchical_cluster_reads(matrix, read_indices=subset_indices)

        # Should return 10 indices (subset of first 10)
        assert len(clustered) == 10
        assert set(clustered) == set(range(10))


class TestSortReads:
    """Tests for the main sort_reads function."""

    def test_sort_reads_no_annotations(self):
        """No annotations should use global clustering."""
        np.random.seed(42)
        matrix = np.random.rand(30, 60)
        names = np.array([f"read_{i:03d}" for i in range(30)])

        indices, sorted_names = sort_reads(matrix, names)

        # Should return all reads
        assert len(indices) == 30
        assert len(sorted_names) == 30
        assert set(sorted_names) == set(names)

        # Indices should be valid
        assert set(indices) == set(range(30))

    def test_sort_reads_single_read(self):
        """Single read should return immediately."""
        matrix = np.array([[0.5, 0.5, 0.5]])
        names = np.array(["read_001"])

        indices, sorted_names = sort_reads(matrix, names)

        assert len(indices) == 1
        assert indices[0] == 0
        assert sorted_names[0] == "read_001"

    def test_sort_reads_qualitative_only(self):
        """Qualitative annotation should group and cluster."""
        np.random.seed(42)
        n_reads = 40
        matrix = np.random.rand(n_reads, 80)
        names = np.array([f"read_{i:03d}" for i in range(n_reads)])

        # Create annotations with haplotype
        haplotypes = np.array(["H1"] * 20 + ["H2"] * 20)
        annot_df = pl.DataFrame({"read_name": names, "haplotype": haplotypes})
        annot = ReadAnnotations(annot_df)

        indices, sorted_names = sort_reads(
            matrix, names, annotations=annot, sort_by=["haplotype"]
        )

        # Should return all reads
        assert len(indices) == n_reads
        assert set(sorted_names) == set(names)

        # Check grouping: H1 reads should be contiguous, then H2
        sorted_haplotypes = [
            annot.data.filter(pl.col("read_name") == name)["haplotype"][0]
            for name in sorted_names
        ]

        # Count transitions between H1 and H2
        transitions = sum(
            1
            for i in range(1, len(sorted_haplotypes))
            if sorted_haplotypes[i] != sorted_haplotypes[i - 1]
        )

        # Should have at most 1 transition (H1 -> H2 or H2 -> H1)
        assert transitions <= 1

    def test_sort_reads_quantitative_only(self):
        """Quantitative annotation should sort in ascending order."""
        n_reads = 30
        matrix = np.random.rand(n_reads, 60)
        names = np.array([f"read_{i:03d}" for i in range(n_reads)])

        # Create annotations with mapq
        np.random.seed(42)
        mapq_values = np.random.randint(10, 60, n_reads)
        annot_df = pl.DataFrame({"read_name": names, "mapq": mapq_values})
        annot = ReadAnnotations(annot_df)

        indices, sorted_names = sort_reads(
            matrix, names, annotations=annot, sort_by=["mapq"]
        )

        # Should return all reads
        assert len(indices) == n_reads

        # Extract mapq values in sorted order
        sorted_mapq = np.array(
            [
                annot.data.filter(pl.col("read_name") == name)["mapq"][0]
                for name in sorted_names
            ]
        )

        # Should be sorted in ascending order
        assert np.all(sorted_mapq[:-1] <= sorted_mapq[1:])

    def test_sort_reads_qualitative_and_quantitative(self):
        """Qualitative + quantitative should group then sort within groups."""
        n_reads = 40
        matrix = np.random.rand(n_reads, 80)
        names = np.array([f"read_{i:03d}" for i in range(n_reads)])

        # Create annotations
        np.random.seed(42)
        haplotypes = np.array(["H1"] * 20 + ["H2"] * 20)
        mapq_values = np.random.randint(10, 60, n_reads)
        annot_df = pl.DataFrame(
            {"read_name": names, "haplotype": haplotypes, "mapq": mapq_values}
        )
        annot = ReadAnnotations(annot_df)

        indices, sorted_names = sort_reads(
            matrix, names, annotations=annot, sort_by=["haplotype", "mapq"]
        )

        # Should return all reads
        assert len(indices) == n_reads

        # Extract haplotypes and mapq in sorted order
        sorted_data = [
            (
                annot.data.filter(pl.col("read_name") == name)["haplotype"][0],
                annot.data.filter(pl.col("read_name") == name)["mapq"][0],
            )
            for name in sorted_names
        ]
        sorted_haplotypes = [h for h, _ in sorted_data]
        sorted_mapq = [m for _, m in sorted_data]

        # Check grouping by haplotype
        transitions = sum(
            1
            for i in range(1, len(sorted_haplotypes))
            if sorted_haplotypes[i] != sorted_haplotypes[i - 1]
        )
        assert transitions <= 1

        # Within each haplotype group, mapq should be sorted
        h1_mapq = [m for h, m in sorted_data if h == "H1"]
        h2_mapq = [m for h, m in sorted_data if h == "H2"]

        assert all(h1_mapq[i] <= h1_mapq[i + 1] for i in range(len(h1_mapq) - 1))
        assert all(h2_mapq[i] <= h2_mapq[i + 1] for i in range(len(h2_mapq) - 1))

    def test_sort_reads_unknown_column(self):
        """Unknown column in sort_by should raise error."""
        matrix = np.random.rand(10, 20)
        names = np.array([f"read_{i:03d}" for i in range(10)])

        annot_df = pl.DataFrame(
            {"read_name": names, "haplotype": ["H1"] * 5 + ["H2"] * 5}
        )
        annot = ReadAnnotations(annot_df)

        with pytest.raises(ValueError, match="unknown annotation columns"):
            sort_reads(matrix, names, annotations=annot, sort_by=["nonexistent"])

    def test_sort_reads_missing_annotations(self):
        """Reads without annotations should raise error."""
        matrix = np.random.rand(20, 40)
        names = np.array([f"read_{i:03d}" for i in range(20)])

        # Only annotate first 10 reads
        annot_df = pl.DataFrame(
            {"read_name": names[:10], "haplotype": ["H1"] * 5 + ["H2"] * 5}
        )
        annot = ReadAnnotations(annot_df)

        with pytest.raises(ValueError, match="have no annotations"):
            sort_reads(matrix, names, annotations=annot, sort_by=["haplotype"])

    def test_sort_reads_empty_sort_by(self):
        """Empty sort_by should use clustering."""
        np.random.seed(42)
        matrix = np.random.rand(15, 30)
        names = np.array([f"read_{i:03d}" for i in range(15)])

        annot_df = pl.DataFrame(
            {"read_name": names, "haplotype": ["H1"] * 8 + ["H2"] * 7}
        )
        annot = ReadAnnotations(annot_df)

        # Empty sort_by list
        indices, sorted_names = sort_reads(
            matrix, names, annotations=annot, sort_by=[]
        )

        # Should return all reads (using clustering, not annotation sorting)
        assert len(indices) == 15
        assert set(sorted_names) == set(names)

    def test_sort_reads_preserves_matrix_alignment(self):
        """Returned indices should correctly align with matrix rows."""
        np.random.seed(42)
        n_reads = 25
        matrix = np.random.rand(n_reads, 50)
        names = np.array([f"read_{i:03d}" for i in range(n_reads)])

        # Add distinct patterns to verify alignment
        matrix[0, :] = 0.0  # First read all unmethylated
        matrix[-1, :] = 1.0  # Last read all methylated

        annot_df = pl.DataFrame(
            {"read_name": names, "haplotype": ["H1"] * 13 + ["H2"] * 12}
        )
        annot = ReadAnnotations(annot_df)

        indices, sorted_names = sort_reads(
            matrix, names, annotations=annot, sort_by=["haplotype"]
        )

        # Verify that sorted_names[i] corresponds to matrix[indices[i]]
        for i, sorted_name in enumerate(sorted_names):
            original_idx = indices[i]
            assert names[original_idx] == sorted_name

            # Verify matrix alignment
            if sorted_name == "read_000":
                assert np.allclose(matrix[original_idx], 0.0)
            elif sorted_name == f"read_{n_reads-1:03d}":
                assert np.allclose(matrix[original_idx], 1.0)
