from typing import List


def cross_check_genes(
    query_genes: List[str], reference_genes: List[str], min_gene_overlap: int = 100
) -> List[str]:
    """Cross check the overlap between two lists of genes

    Parameters
    ----------
    query_genes : List[str]
        List of gene names
    reference_genes : List[str]
       List of gene names

    Returns
    -------
    List[str]
        List of overlapping genes
    """
    common_ens_ids = list(set(reference_genes).intersection(set(query_genes)))
    assert len(common_ens_ids) >= min_gene_overlap, (
        f"The intersection of genes between the query and reference dataset is too small, expected at least {min_gene_overlap}."
    )

    return common_ens_ids
