#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import warnings

import numpy as np
import scipy.sparse as sp
import scipy.stats as st


def library_size(
    n_cells: int,
    loc: float = 7.5,
    scale: float = 0.5,
    lower_bound: float = -1.0,
    upper_bound: float = np.inf,
) -> np.ndarray:
    """log-normal noise for the number of umis per cell (with a lower bound to
    represent a minimum depth cutoff)

    :param n_cells: number of library sizes to generate
    :param loc: mean of library size in log-space
    :param scale: standard deviation
    :param lower_bound: lower bound relative to ``loc``
    :param upper_bound: upper bound relative to ``loc``
    :return: float array of shape (n_cells,) containing the library sizes
    """

    return np.exp(
        st.truncnorm.rvs(lower_bound, upper_bound, loc=loc, scale=scale, size=n_cells)
    ).astype(int)


def fragment_genes(n_genes: int, lam: float = 1.0) -> np.ndarray:
    """Generate a random number of fragments for each gene used a poisson distribution

    :param n_genes: number of genes to fragment
    :param lam: mean of poisson distribution of (additional) fragments. Every gene
                will have at least one fragment
    :return: int array of shape (n_genes,) with number of fragments per gene
    """
    # random number of possible fragments per gene, poisson distributed
    # add one to ensure ≥1 fragment per gene
    fragments_per_gene = 1 + np.random.poisson(lam, size=n_genes)

    return fragments_per_gene


def umi_counts(
    raw_expression: np.ndarray,
    lib_size: int | np.ndarray[int] = None,
    fragments_per_gene: int | np.ndarray[int] = 1,
    sparse: bool = False,
) -> np.ndarray:
    """Given an ``(..., n_genes)`` array of expression values, generates a count matrix
    of UMIs based by multinomial sampling. The last dimension of the array defines
    the number of genes per cell, while all other dimensions are assumed to index
    distinct cells.

    :param raw_expression: array of raw expression values (non-negative)
    :param lib_size: library size for each cell, either constant or per-sample.
                     If None, generates a distribution using ``library_size``
    :param fragments_per_gene: fragments observed per gene, either constant or per-gene
    :param sparse: if True, return a sparse (CSR) matrix. If ``raw_expression`` is >2-d,
                   the output matrix will be reshaped to 2-d
    :return: integer array of shape ``(..., n_features)`` containing umi counts
    """

    if np.any(raw_expression < 0):
        raise ValueError("raw_expression must be non-negative")

    if len(raw_expression.shape) == 1:
        raise ValueError("raw_expression should be >= 2 dimensions")

    if sparse and len(raw_expression.shape) > 2:
        warnings.warn("Reshaping output to 2 dimensions for sparse matrix")

    n_cells = raw_expression.shape[:-1]
    n_genes = raw_expression.shape[-1]

    if lib_size is None:
        lib_size = library_size(n_cells)  # is this a good default?
    else:
        lib_size = np.broadcast_to(lib_size, n_cells)

    fragments_per_gene = np.broadcast_to(fragments_per_gene, (n_genes,))

    # each fragment is at the level of the gene it comes from
    fragment_expression = np.repeat(raw_expression, fragments_per_gene, axis=-1)

    gene_p = fragment_expression / fragment_expression.sum(-1, keepdims=True)

    if sparse:
        cell_gene_umis = sp.vstack(
            [
                sp.csr_matrix(np.random.multinomial(n=lib_size[i], pvals=gene_p[i]))
                for i in np.ndindex(*n_cells)
            ]
        )
    else:
        cell_gene_umis = np.vstack(
            [
                np.random.multinomial(n=lib_size[i], pvals=gene_p[i])
                for i in np.ndindex(*n_cells)
            ]
        ).reshape(*n_cells, -1)

    return cell_gene_umis


def pcr_noise(
    read_counts: np.ndarray,
    pcr_betas: float | np.ndarray[float],
    n_cycles: int,
    copy: bool = True,
) -> np.ndarray:
    """PCR noise model: every read has an affinity for PCR, and for every round
    of PCR we do a ~binomial doubling of each count.

    :param read_counts: array of shape (n_samples, n_features) representing unique
                        molecules (e.g. genes or gene fragments). If a sparse matrix
                        is provided, the output will be sparse
    :param pcr_betas: PCR efficiency for each feature, either constant or per-feature
    :param n_cycles: number of rounds of PCR to simulate
    :param copy: if True, return a copy of the read_counts array, else modify in-place
    :return: int array of shape (n_samples, n_features) with amplified counts
    """
    if np.any(pcr_betas < 0):
        raise ValueError("pcr_betas must be non-negative")

    if copy:
        read_counts = read_counts.copy()

    pcr_betas = np.broadcast_to(pcr_betas, (1, read_counts.shape[1]))

    if sp.issparse(read_counts):
        d = read_counts.data[None, :]
        pcr_betas = pcr_betas[:, read_counts.nonzero()[1]]
    else:
        d = read_counts

    # for each round of pcr, each gene increases according to its affinity factor
    for i in range(n_cycles):
        d += np.random.binomial(n=d, p=pcr_betas, size=d.shape)

    if sp.issparse(read_counts):
        read_counts.data = d.flatten()

    return read_counts
