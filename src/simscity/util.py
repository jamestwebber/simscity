#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import warnings

import numpy as np
import sparse


def arrays_to_anndata(
    expression: np.ndarray | sparse.GCXS,
    batch: np.ndarray,
    classes=np.ndarray,
    **obsm: np.ndarray,
):
    """Packages together an AnnData object for use in existing pipelines

    :param expression: array of gene expression
    :param batch: array of batch assignments
    :param classes: array of class (e.g. cell type) assignments
    :param obsm: any other arrays to store as metadata
    :return: AnnData object containing the given data and metadata
    """
    try:
        import anndata
        import pandas as pd
    except ImportError:
        warnings.warn("arrays_to_anndata requires anndata")
        raise

    metadata = pd.DataFrame({"batch": batch, "class": classes})
    metadata["batch"] = metadata["batch"].astype("category")
    metadata["class"] = metadata["class"].astype("category")

    adata = anndata.AnnData(X=expression, obs=metadata, obsm=(obsm or None))

    return adata
