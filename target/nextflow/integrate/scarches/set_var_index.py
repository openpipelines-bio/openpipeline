import anndata as ad


def strip_version_number(gene_series: list[str]) -> list[str]:
    """Sanitize ensemble ID's by removing version numbers.

    Parameters
    ----------
    gene_series : list[str]
        List of ensemble ID's to sanitize.

    Returns
    -------
    list[str]
        List of sanitized ensemble ID's.
    """

    # Convert to string type to handle Categorical series
    gene_series = gene_series.astype(str)

    # Pattern matches Ensembl IDs: starts with ENS, followed by any characters,
    # then an eleven digit number, optionally followed by .version_number
    ensembl_pattern = r"^(ENS.*\d{11})(?:\.\d+)?$"
    ensembl_mask = gene_series.str.match(ensembl_pattern)

    sanitized = gene_series.where(
        ~ensembl_mask, gene_series.str.extract(ensembl_pattern)[0]
    )

    return sanitized


def set_var_index(
    adata: ad.AnnData, var_name: str | None = None, sanitize_ensembl_ids: bool = True
) -> ad.AnnData:
    """Sanitize gene names (optional) and set the index of the .var DataFrame.

    Parameters
    ----------
    adata : AnnData
        Annotated data object
    var_name : str | None
        Name of the column in `adata.var` that contains the gene names, if None, the existing index will be sanitized but not replaced.
    sanitize_ensembl_ids : bool
        Whether to sanitize gene names by removing version numbers.

    Returns
    -------
    AnnData
        Copy of `adata` with optionally sanitized and replaced index
    """
    gene_names = adata.var[var_name] if var_name else adata.var.index.to_series()

    if sanitize_ensembl_ids:
        ori_gene_names = len(gene_names)
        gene_names = strip_version_number(gene_names)
        sanitized_gene_names = len(set(gene_names))

        assert ori_gene_names == sanitized_gene_names, (
            "Sanitizing gene names resulted in duplicated gene names.\n"
            "Please ensure unique gene names before proceeding.\n"
            "Please make sure --var_gene_names contains ensembl IDs (not gene symbols) "
            "when --sanitize_ensembl_ids is set to True."
        )

    adata.var.index = gene_names

    return adata
