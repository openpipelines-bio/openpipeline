def subset_vars(adata, subset_col):
    """Subset AnnData object on highly variable genes
    
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
    if not subset_col in adata.var.columns:
        raise ValueError(f"Requested to use .var column '{subset_col}' as a selection of genes, but the column is not available.")

    return adata[:, adata.var[subset_col]].copy()
