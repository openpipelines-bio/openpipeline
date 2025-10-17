import anndata as ad
import re


def strip_version_number(gene_names: list[str]) -> list[str]:
    """Sanitize gene names by removing version numbers.

    Parameters
    ----------
    gene_names : list[str]
        List of gene names to sanitize.

    Returns
    -------
    list[str]
        List of sanitized gene names.
    """
    return [re.sub("\\.[0-9]+$", "", s) for s in gene_names]


def set_var_index(
    adata: ad.AnnData, var_name: str | None = None, sanitize_gene_names: bool = True
) -> ad.AnnData:
    """Sanitize gene names (optional) and set the index of the .var DataFrame.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    var_name : str | None
        Name of the column in `adata.var` that contains the gene names, if None, the existing index will be sanitized but not replaced.
    sanitize_gene_names : bool
        Whether to sanitize gene names by removing version numbers.

    Returns
    -------
    AnnData
        Copy of `adata` with optionally sanitized and replaced index
    """
    gene_names = adata.var[var_name] if var_name else adata.var.index

    if sanitize_gene_names:
        ori_gene_names = len(gene_names)
        gene_names = strip_version_number(gene_names)
        sanitized_gene_names = len(set(gene_names))

        assert ori_gene_names == sanitized_gene_names, (
            "Sanitizing gene names resulted in duplicated gene names.\n"
            "Please ensure unique gene names before proceeding.\n"
            "Please make sure --var_gene_names contains ensembl IDs (not gene symbols) "
            "when --sanitize_gene_names is set to True."
        )

    adata.var.index = gene_names

    return adata
