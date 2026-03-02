"""
Read sorting and clustering algorithms.

This module implements read ordering strategies including hierarchical clustering
based on methylation patterns and multi-column annotation-based sorting.
"""

from typing import Optional
import numba
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list

from cpgplotter.core.annotations import ReadAnnotations


def _sort_by_clustering(
    methylation_matrix: np.ndarray,
    read_names: np.ndarray,
    nan_weight: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sort reads using global hierarchical clustering.

    Args:
        methylation_matrix: Methylation probability matrix (n_reads x n_cpgs)
        read_names: Array of read names
        nan_weight: Maximum weight for coverage distance component (0-1)

    Returns:
        Tuple of (sorted_indices, sorted_read_names)
    """
    indices = hierarchical_cluster_reads(methylation_matrix, nan_weight=nan_weight)
    return indices, read_names[indices]


def _sort_by_qualitative_and_cluster(
    methylation_matrix: np.ndarray,
    read_names: np.ndarray,
    annotations: ReadAnnotations,
    qualitative_cols: list[str],
    nan_weight: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Group reads by qualitative annotations and cluster within each group.

    Args:
        methylation_matrix: Methylation probability matrix (n_reads x n_cpgs)
        read_names: Array of read names
        annotations: ReadAnnotations object
        qualitative_cols: List of qualitative column names to group by
        nan_weight: Maximum weight for coverage distance component (0-1)

    Returns:
        Tuple of (sorted_indices, sorted_read_names)
    """
    # Create read_name -> matrix index mapping
    read_name_to_idx = {name: idx for idx, name in enumerate(read_names)}

    result_indices = []

    # Sort annotations by qualitative columns for deterministic group order
    sorted_annot = annotations.data.sort(qualitative_cols)

    # Group by qualitative columns
    for group_key, group_df in sorted_annot.group_by(
        qualitative_cols, maintain_order=True
    ):
        # Get read names in this group
        group_read_names = group_df["read_name"].to_numpy()

        # Map to matrix indices
        group_indices = np.array([read_name_to_idx[name] for name in group_read_names])

        # Apply hierarchical clustering within group
        clustered_indices = hierarchical_cluster_reads(
            methylation_matrix, read_indices=group_indices, nan_weight=nan_weight
        )

        result_indices.extend(clustered_indices)

    result_indices = np.array(result_indices)
    return result_indices, read_names[result_indices]


def _sort_by_qualitative_and_quantitative(
    read_names: np.ndarray,
    annotations: ReadAnnotations,
    qualitative_cols: list[str],
    sort_by: list[str],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Group reads by qualitative annotations and sort within groups by quantitative.

    Args:
        read_names: Array of read names
        annotations: ReadAnnotations object
        qualitative_cols: List of qualitative column names to group by
        sort_by: List of all columns to sort by (qualitative + quantitative)

    Returns:
        Tuple of (sorted_indices, sorted_read_names)
    """
    # Create read_name -> matrix index mapping
    read_name_to_idx = {name: idx for idx, name in enumerate(read_names)}

    result_indices = []

    # Sort annotations by all sort_by columns
    sorted_annot = annotations.data.sort(sort_by)

    # Group by qualitative columns
    for group_key, group_df in sorted_annot.group_by(
        qualitative_cols, maintain_order=True
    ):
        # Within each group, reads are already sorted by quantitative cols
        group_read_names = group_df["read_name"].to_numpy()

        # Map to matrix indices
        group_indices = np.array([read_name_to_idx[name] for name in group_read_names])

        result_indices.extend(group_indices)

    result_indices = np.array(result_indices)
    return result_indices, read_names[result_indices]


def _sort_by_quantitative(
    read_names: np.ndarray, annotations: ReadAnnotations, quantitative_cols: list[str]
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sort reads by quantitative annotations only.

    Args:
        read_names: Array of read names
        annotations: ReadAnnotations object
        quantitative_cols: List of quantitative column names to sort by

    Returns:
        Tuple of (sorted_indices, sorted_read_names)
    """
    # Create read_name -> matrix index mapping
    read_name_to_idx = {name: idx for idx, name in enumerate(read_names)}

    # Sort entire dataset by quantitative columns
    sorted_annot = annotations.data.sort(quantitative_cols)
    sorted_read_names = sorted_annot["read_name"].to_numpy()

    # Map to matrix indices
    sorted_indices = np.array([read_name_to_idx[name] for name in sorted_read_names])

    return sorted_indices, read_names[sorted_indices]


def sort_reads(
    methylation_matrix: np.ndarray,
    read_names: np.ndarray,
    annotations: Optional[ReadAnnotations] = None,
    sort_by: Optional[list[str]] = None,
    nan_weight: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Sort reads according to specified criteria.

    Sorting priority:
    1. If sort_by is specified: sort by those annotation columns in order
    2. If no quantitative annotations: use hierarchical clustering within groups
    3. Default: hierarchical clustering on all reads

    Args:
        methylation_matrix: Methylation probability matrix (n_reads x n_cpgs)
        read_names: Array of read names
        annotations: Optional ReadAnnotations object
        sort_by: Optional list of annotation column names for sorting
        nan_weight: Maximum weight for coverage distance component in
            clustering (0-1). At 0 the metric is pure Hellinger distance.
            Default: 0.5.

    Returns:
        Tuple of (sorted_indices, sorted_read_names)
        - sorted_indices: Array of indices for reordering
        - sorted_read_names: Reordered read names
    """
    n_reads = len(read_names)

    # Handle edge cases
    if n_reads == 0:
        return np.array([], dtype=int), read_names

    if n_reads == 1:
        return np.array([0]), read_names

    # Case 1: No annotations or no sort criteria - use global clustering
    if annotations is None or sort_by is None or len(sort_by) == 0:
        return _sort_by_clustering(methylation_matrix, read_names, nan_weight=nan_weight)

    # Validate sort_by columns exist
    missing_cols = set(sort_by) - set(annotations.column_types.keys())
    if missing_cols:
        available_cols = list(annotations.column_types.keys())
        raise ValueError(
            f"sort_by contains unknown annotation columns: {missing_cols}\n"
            f"Available columns: {available_cols}"
        )

    # Filter annotations to reads in matrix
    aligned_annot = annotations.filter_reads(set(read_names))

    # Check all matrix reads have annotations
    annot_names = set(aligned_annot.data["read_name"])
    missing = set(read_names) - annot_names
    if missing:
        example_missing = list(missing)[:5]
        raise ValueError(
            f"Cannot sort: {len(missing)} reads in methylation matrix "
            f"have no annotations. All reads must be annotated when using sort_by.\n"
            f"Example missing reads: {example_missing}"
        )

    # Separate qualitative and quantitative columns in sort_by
    qualitative_cols = [
        col for col in sort_by if annotations.column_types[col] == "qualitative"
    ]
    quantitative_cols = [
        col for col in sort_by if annotations.column_types[col] == "quantitative"
    ]

    # Route to appropriate sorting function based on column types
    if len(quantitative_cols) == 0:
        # Only qualitative columns (or no valid columns)
        if len(qualitative_cols) == 0:
            return _sort_by_clustering(methylation_matrix, read_names, nan_weight=nan_weight)
        return _sort_by_qualitative_and_cluster(
            methylation_matrix, read_names, aligned_annot, qualitative_cols,
            nan_weight=nan_weight,
        )
    elif len(qualitative_cols) > 0:
        # Both qualitative and quantitative columns
        return _sort_by_qualitative_and_quantitative(
            read_names, aligned_annot, qualitative_cols, sort_by
        )
    else:
        # Only quantitative columns
        return _sort_by_quantitative(read_names, aligned_annot, quantitative_cols)


def hierarchical_cluster_reads(
    methylation_matrix: np.ndarray,
    read_indices: Optional[np.ndarray] = None,
    nan_weight: float = 0.5,
) -> np.ndarray:
    """
    Cluster reads by methylation pattern using hierarchical clustering.

    Uses a coverage-aware Hellinger distance that combines methylation
    similarity with coverage pattern similarity. The coverage component
    is weighted by a sigmoid of the pair NaN fraction.

    Args:
        methylation_matrix: Methylation probability matrix (n_reads x n_cpgs)
        read_indices: Optional subset of read indices to cluster.
            If None, clusters all reads.
        nan_weight: Maximum weight for coverage distance component (0-1).

    Returns:
        Array of indices representing the clustering order
    """
    if read_indices is None:
        read_indices = np.arange(len(methylation_matrix))

    subset = methylation_matrix[read_indices]

    # Compute coverage-aware distance for reads
    distances = compute_distance_matrix(subset, nan_weight=nan_weight)

    # Perform hierarchical clustering with Ward's method
    linkage_matrix = linkage(distances, method="ward")

    # Get leaf order
    order = leaves_list(linkage_matrix)

    return read_indices[order]


def compute_distance_matrix(
    methylation_matrix: np.ndarray,
    nan_weight: float = 0.5,
) -> np.ndarray:
    """
    Compute pairwise distances between reads using a coverage-aware Hellinger metric.

    Combines two components:
    1. Normalized Hellinger distance over CpGs covered by both reads (methylation similarity)
    2. Jaccard distance on coverage patterns (NaN pattern similarity)

    The contribution of the coverage component is modulated by a sigmoid function
    of the pair's NaN fraction, so it has negligible effect when reads are well-covered
    (<10% NaN) and dominates when reads are sparse (>30% NaN).

    Args:
        methylation_matrix: Methylation probability matrix (n_reads x n_cpgs)
        nan_weight: Maximum weight for the coverage distance component (0-1).
            At 0 the metric is pure Hellinger. Default: 0.5.

    Returns:
        Condensed distance matrix (1D array of length n*(n-1)/2)
    """
    return _compute_distances_numba(
        np.ascontiguousarray(methylation_matrix, dtype=np.float64),
        nan_weight,
    )


# Keep old name as alias for backwards compatibility in tests
compute_hellinger_distance_matrix = compute_distance_matrix


@numba.njit
def _sigmoid(x: float, k: float, midpoint: float) -> float:
    """Logistic sigmoid: 1 / (1 + exp(-k * (x - midpoint)))."""
    return 1.0 / (1.0 + np.exp(-k * (x - midpoint)))


@numba.njit
def _compute_distances_numba(
    matrix: np.ndarray, nan_weight: float
) -> np.ndarray:
    """Numba-accelerated pairwise coverage-aware Hellinger distance."""
    n_reads = matrix.shape[0]
    n_cpgs = matrix.shape[1]
    no_overlap_distance = 1.0  # max normalised distance

    # Sigmoid parameters for non-linear NaN weighting
    sigmoid_k = 20.0
    sigmoid_mid = 0.2

    # Pre-compute per-read NaN fractions
    nan_fracs = np.empty(n_reads, dtype=np.float64)
    for i in range(n_reads):
        n_nan = 0
        for k in range(n_cpgs):
            if np.isnan(matrix[i, k]):
                n_nan += 1
        nan_fracs[i] = n_nan / n_cpgs if n_cpgs > 0 else 0.0

    condensed_size = n_reads * (n_reads - 1) // 2
    distances = np.empty(condensed_size, dtype=np.float64)

    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    idx = 0
    for i in range(n_reads):
        for j in range(i + 1, n_reads):
            sum_sq = 0.0
            n_overlap = 0  # both covered
            n_both_covered_or_nan = 0  # |C_i ∩ C_j| (both not-NaN)
            n_either_covered = 0  # |C_i ∪ C_j| (at least one not-NaN)

            for k in range(n_cpgs):
                vi = matrix[i, k]
                vj = matrix[j, k]
                i_nan = np.isnan(vi)
                j_nan = np.isnan(vj)

                if not i_nan and not j_nan:
                    # Both covered — contributes to Hellinger
                    diff = np.sqrt(vi) - np.sqrt(vj)
                    sum_sq += diff * diff
                    n_overlap += 1
                    n_both_covered_or_nan += 1
                    n_either_covered += 1
                elif i_nan and j_nan:
                    # Both NaN — agree on coverage (not in union)
                    n_both_covered_or_nan += 1
                else:
                    # One covered, one NaN — disagree on coverage
                    n_either_covered += 1

            # Component 1: normalised Hellinger distance [0, 1]
            if n_overlap > 0:
                hellinger = np.sqrt(0.5 * sum_sq / n_overlap) * inv_sqrt2
                # Clamp to [0, 1] for numerical safety
                if hellinger > 1.0:
                    hellinger = 1.0
            else:
                hellinger = no_overlap_distance

            # Component 2: Jaccard distance on coverage [0, 1]
            if n_either_covered > 0:
                jaccard = 1.0 - n_overlap / n_either_covered
            else:
                # Both reads entirely NaN — identical coverage pattern
                jaccard = 0.0

            # Non-linear weight: sigmoid of max(nan_frac_i, nan_frac_j)
            pair_nan_frac = max(nan_fracs[i], nan_fracs[j])
            w = nan_weight * _sigmoid(pair_nan_frac, sigmoid_k, sigmoid_mid)

            distances[idx] = (1.0 - w) * hellinger + w * jaccard
            idx += 1

    return distances
