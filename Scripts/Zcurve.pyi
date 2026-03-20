"""
Zcurve module for DNA sequence feature extraction using Z-curve method.

This module provides functionality to encode DNA sequences into numerical feature
vectors using the Z-curve method, which is a powerful technique for representing
genomic sequences in a three-dimensional space. The Z-curve method captures
the compositional properties of DNA sequences and is particularly useful for
gene prediction and sequence analysis.

Functions:
    encode(records, n_jobs=0) - Encode DNA sequences into Z-curve feature vectors

Example:
    >>> import Zcurve
    >>> sequences = ["ATCGATCG", "GCTAGCTA"]
    >>> features = Zcurve.encode(sequences)
    >>> print(features.shape)
    (2, 765)

The module is optimized for performance using AVX/FMA instructions and supports
parallel processing through multi-threading for efficient handling of large
datasets.
"""
import numpy as np
from typing import List

def encode(records: List[str], n_jobs: int = 0) -> np.ndarray:
    """
    Encode DNA sequences into Z-curve feature vectors.
    
    Parameters:
        records (List[str]): List of DNA sequences as strings. Each string should
            contain valid DNA nucleotides (A, T, C, G and other degenerate bases).
        n_jobs (int, optional): Number of threads to use for parallel processing.
            Default is 0, which means using all available cores.
    
    Returns:
        np.ndarray: A 2D array where each row represents the Z-curve features
            for a sequence. The number of columns depends on the feature extraction
            method (typically 765 features per sequence).
    
    Raises:
        TypeError: If records is not a list of strings.
    
    Example:
        >>> import Zcurve
        >>> # Encode a single sequence
        >>> features = Zcurve.encode(["ATCGATCGATCG"])
        >>> # Encode multiple sequences with parallel processing
        >>> features = Zcurve.encode(["ATCG", "GCTA"], n_jobs=2)
    """
    ...