import anndata as ad
import re


def set_var_index(adata: ad.AnnData, var_name: str | None = None) -> ad.AnnData:
    """Sanitize gene names and set the index of the .var DataFrame.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    var_name : str | None
        Name of the column in `adata.var` that contains the gene names, if None, the existing index will be sanitized but not replaced.

    Returns
    -------
    AnnData
        Copy of `adata` with sanitized and replaced index
    """
    if var_name:
        adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in adata.var[var_name]]
    else:
        adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in adata.var.index]
    return adata
