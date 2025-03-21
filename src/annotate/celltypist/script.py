import sys
import celltypist
import mudata as mu
import anndata as ad
import pandas as pd
import numpy as np

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
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
    "input_var_gene_names": None,
    "reference_var_gene_names": "ensemblid",
    "input_layer": "log_normalized",
    "reference_layer": None,
    "output_obs_predictions": "celltypist_pred",
    "output_obs_probability": "celltypist_probability",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from cross_check_genes import cross_check_genes
from set_var_index import set_var_index
from subset_vars import subset_vars

logger = setup_logger()


def check_celltypist_format(indata):
    if np.abs(np.expm1(indata[0]).sum() - 10000) > 1:
        return False
    return True


def main(par):
    if (not par["model"] and not par["reference"]) or (
        par["model"] and par["reference"]
    ):
        raise ValueError(
            "Make sure to provide either 'model' or 'reference', but not both."
        )

    input_mudata = mu.read_h5mu(par["input"])
    input_adata = input_mudata.mod[par["modality"]]
    input_modality = input_adata.copy()

    # Provide correct format of query data for celltypist annotation
    ## Sanitize gene names and set as index
    input_modality = set_var_index(input_modality, par["input_var_gene_names"])
    ## Fetch lognormalized counts
    lognorm_counts = (
        input_modality.layers[par["input_layer"]].copy()
        if par["input_layer"]
        else input_modality.X.copy()
    )
    ## Create AnnData object
    input_modality = ad.AnnData(
        X=lognorm_counts, var=pd.DataFrame(index=input_modality.var.index)
    )

    if par["model"]:
        logger.info("Loading CellTypist model")
        model = celltypist.models.Model.load(par["model"])
        cross_check_genes(
            input_modality.var.index,
            model.features,
            min_gene_overlap=par["input_reference_gene_overlap"],
        )

    elif par["reference"]:
        reference_modality = mu.read_h5mu(par["reference"]).mod[par["modality"]]

        # subset to HVG if required
        if par["reference_var_input"]:
            reference_modality = subset_vars(
                reference_modality, par["reference_var_input"]
            )

        # Set var names to the desired gene name format (gene symbol, ensembl id, etc.)
        # CellTypist requires query gene names to be in index
        reference_modality = set_var_index(
            reference_modality, par["reference_var_gene_names"]
        )

        # Ensure enough overlap between genes in query and reference
        cross_check_genes(
            input_modality.var.index,
            reference_modality.var.index,
            min_gene_overlap=par["input_reference_gene_overlap"],
        )

        reference_matrix = (
            reference_modality.layers[par["reference_layer"]]
            if par["reference_layer"]
            else reference_modality.X
        )

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
            check_expression=True,
        )

    logger.info("Predicting CellTypist annotations")
    predictions = celltypist.annotate(
        input_modality, model, majority_voting=par["majority_voting"]
    )
    input_adata.obs[par["output_obs_predictions"]] = predictions.predicted_labels[
        "predicted_labels"
    ]
    input_adata.obs[par["output_obs_probability"]] = predictions.probability_matrix.max(
        axis=1
    ).values

    # copy observations back to input data (with full set of features)
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main(par)
