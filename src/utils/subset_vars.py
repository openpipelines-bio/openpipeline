def subset_vars(adata, subset_col):
    """
    Subset AnnData object on highly variable genes or a boolean mask.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    subset_col : str, pd.Series, pd.Index, or np.ndarray
        Name of the boolean column in `adata.var` that contains the information if features should be used or not,
        or a boolean mask (same length as adata.var)

    Returns
    -------
    AnnData
        Copy of `adata` with subsetted features
    """
    import pandas as pd
    import numpy as np

    if isinstance(subset_col, str):
        if subset_col not in adata.var.columns:
            raise ValueError(
                f"Requested to use .var column '{subset_col}' as a selection of genes, but the column is not available."
            )
        mask = adata.var[subset_col]
        if mask.dtype == "boolean":
            assert mask.isna().sum() == 0, (
                f"The .var column `{subset_col}` contains NaN values. Can not subset data."
            )
            mask = mask.astype("bool")
        assert mask.dtype == "bool", (
            f"Expected dtype of .var column '{subset_col}' to be `bool`, but found {mask.dtype}. Can not subset data."
        )
    elif isinstance(subset_col, (pd.Series, pd.Index, np.ndarray, list)):
        mask = pd.Series(subset_col, index=adata.var.index) if not isinstance(subset_col, pd.Series) else subset_col
        if mask.dtype != "bool":
            raise ValueError(f"Provided mask must be boolean, got {mask.dtype}.")
        if mask.isna().sum() > 0:
            raise ValueError("Provided mask contains NaN values. Can not subset data.")
        if len(mask) != adata.n_vars:
            raise ValueError(f"Mask length {len(mask)} does not match number of variables {adata.n_vars}.")
    else:
        raise TypeError("subset_col must be a string (column name) or a boolean mask (Series, Index, ndarray, or list).")

    return adata[:, mask].copy()
