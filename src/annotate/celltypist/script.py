import sys
import logging
import celltypist
import mudata as mu
import re
import numpy as np

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix_log_normalized.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": None,
    # "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "model": "resources_test/annotation_test_data/celltypist_model_Immune_All_Low.pkl",
    "input_reference_gene_overlap": 100,
    "reference_obs_target": "cell_ontology_class",
    "reference_var_input": None,
    "check_expression": False,
    "feature_selection": True,
    "majority_voting": True,
    "output_compression": "gzip",
    "input_var_gene_names": "gene_symbol",
    "reference_var_gene_names": "ensemblid",
    "input_layer": None,
    "reference_layer": None,
    "output_obs_predictions": "celltypist_pred",
    "output_obs_probabilities": "celltypist_probability",
}
meta = {
}
## VIASH END

# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
def setup_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger

# from query_reference_allignment import set_var_index, cross_check_genes, subset_vars
logger = setup_logger()
import anndata as ad
from typing import List

def set_var_index(adata: ad.AnnData, var_name: str | None = None):
    if var_name:
        adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in adata.var[var_name]]
    else:
        adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in adata.var.index]
    return adata

def cross_check_genes(query_genes: List[str], reference_genes: List[str], min_gene_overlap: int = 100):
    logger.info("Detecting common vars based on gene ids")
    common_ens_ids = list(set(reference_genes).intersection(set(query_genes)))

    logger.info("  reference n_vars: %i", len(reference_genes))
    logger.info("  input n_vars: %i", len(query_genes))
    logger.info("  intersect n_vars: %i", len(common_ens_ids))
    assert len(common_ens_ids) >= min_gene_overlap, "The intersection of genes between the query and reference dataset is too small."

    return common_ens_ids

def subset_vars(adata: ad.AnnData, var_column: str | None = None):
    if var_column:
        return adata[:, adata.var[var_column]]
    else:
        return adata
# END TEMPORARY WORKAROUND setup_logger


def check_celltypist_format(indata):
    if np.abs(np.expm1(indata[0]).sum()-10000) > 1:
        return False
    return True


def main(par):

    if (not par["model"] and not par["reference"]) or (par["model"] and par["reference"]):
        raise ValueError("Make sure to provide either 'model' or 'reference', but not both.")

    input_mudata = mu.read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]].copy()

    # Set var names to the desired gene name format (gene symbol, ensembl id, etc.)
    # CellTypist requires query gene names to be in index
    input_modality = set_var_index(input_modality, par["input_var_gene_names"])

    if par["model"]:
        logger.info("Loading CellTypist model")
        model = celltypist.models.Model.load(par["model"])
        cross_check_genes(input_modality.var.index, model.features, min_gene_overlap=par["input_reference_gene_overlap"])

    elif par["reference"]:
        reference_modality = mu.read_h5mu(par["reference"]).mod[par["modality"]]

        # subset to HVG if required
        reference_modality = subset_vars(reference_modality, par["reference_var_input"])

        # Set var names to the desired gene name format (gene symbol, ensembl id, etc.)  
        # CellTypist requires query gene names to be in index
        reference_modality = set_var_index(reference_modality, par["reference_var_gene_names"])

        # Ensure enough overlap between genes in query and reference
        cross_check_genes(input_modality.var.index, reference_modality.var.index, min_gene_overlap=par["input_reference_gene_overlap"])

        input_matrix = input_modality.layers[par["input_layer"]] if par["input_layer"] else input_modality.X
        reference_matrix = reference_modality.layers[par["reference_layer"]] if par["reference_layer"] else reference_modality.X

        if not check_celltypist_format(input_matrix):
            logger.warning("Input data is not in the reccommended format for CellTypist.")
        if not check_celltypist_format(reference_matrix):
            logger.warning("Reference data is not in the reccommended format for CellTypist.")

        labels = reference_modality.obs[par["reference_obs_target"]]

        logger.info("Training CellTypist model on reference")
        model = celltypist.train(
            reference_matrix,
            labels=labels,
            genes=reference_modality.var.index,
            C=par["C"],
            max_iter=par["max_iter"],
            use_SGD=par["use_SGD"],
            feature_selection=par["feature_selection"],
            check_expression=par["check_expression"]
            )

    logger.info("Predicting CellTypist annotations")
    predictions = celltypist.annotate(
        input_modality,
        model,
        majority_voting=par["majority_voting"]
        )
    input_modality.obs[par["output_obs_predictions"]] = predictions.predicted_labels["predicted_labels"]
    input_modality.obs[par["output_obs_probability"]] = predictions.probability_matrix.max(axis=1).values

    # copy observations back to input data (with full set of features)
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == '__main__':
    main(par)
