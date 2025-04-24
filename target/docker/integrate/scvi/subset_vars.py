def subset_vars(adata, subset_col):
    """Subset highly variable genes from AnnData object
    
    Parameters
    ----------
    adata : AnnData
        Annotated data object
    subset_col : str
        Name of the boolean column in `adata.var` that contains the information if features should be used or not

    Returns
    -------
    AnnData
        Copy of `adata` with subsetted features
    """
    return adata[:, adata.var[subset_col]].copy()
