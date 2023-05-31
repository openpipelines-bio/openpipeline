def subset_hvg(adata, hvg_col):
    """Subset highly variable genes from AnnData object
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object
    hvg_col : str
        Name of the boolean column in `adata.var` that contains the information if genes should be used or not

    Returns
    -------
    AnnData
        Copy of `adata` with subsetted genes
    """
    return adata[:, adata.var[hvg_col]].copy()
