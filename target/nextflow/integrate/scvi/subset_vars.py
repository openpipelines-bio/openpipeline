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
    if subset_col not in adata.var.columns:
        raise ValueError(
            f"Requested to use .var column '{subset_col}' as a selection of genes, but the column is not available."
        )

    if adata.var[subset_col].dtype == "boolean":
        assert adata.var[subset_col].isna().sum() == 0, (
            f"The .var column `{subset_col}` contains NaN values. Can not subset data."
        )
        adata.var[subset_col] = adata.var[subset_col].astype("bool")

    assert adata.var[subset_col].dtype == "bool", (
        f"Expected dtype of .var column '{subset_col}' to be `bool`, but found {adata.var[subset_col].dtype}. Can not subset data."
    )

    return adata[:, adata.var[subset_col]].copy()
