import re

import anndata as ad


def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger


# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()


# Helper functions
def set_var_index(adata: ad.AnnData, var_name: str | None = None):
    if var_name:
        adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in adata.var[var_name]]
    else:
        adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in adata.var.index]
    return adata


def cross_check_genes(query: ad.AnnData, reference: ad.AnnData):
    logger.info("Detecting common vars based on gene ids")
    common_ens_ids = list(set(reference.var.index).intersection(set(query.var.index)))

    logger.info("  reference n_vars: %i", reference.n_vars)
    logger.info("  input n_vars: %i", query.n_vars)
    logger.info("  intersect n_vars: %i", len(common_ens_ids))
    assert (
        len(common_ens_ids) >= 100
    ), "The intersection of genes between the query and reference dataset is too small."

    return common_ens_ids


def subset_vars(adata: ad.AnnData, var_column: str | None = None):
    if var_column:
        return adata[:, adata.var[var_column]]
    else:
        return adata
