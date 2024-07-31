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
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "model": None,
    "reference_obs_targets": "cell_ontology_class",
    "check_expression": False,
    "feature_selection": True,
    "majority_voting": True,
    "output_compression": "gzip",
    "var_query_gene_names": None,
    "var_reference_gene_names": "ensemblid",
    "reference_layer": None,
    "output_obs_predictions": "celltypist_pred",
    "output_obs_probability": "celltypist_probability",
}
meta = {
}
## VIASH END

# START TEMPORARY WORKAROUND setup_logger
def setup_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger

def check_celltypist_format(indata):
    if np.abs(np.expm1(indata[0]).sum()-10000) > 1:
        return False
    return True

def main(par):
    logger = setup_logger()

    input_mudata = mu.read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]].copy()
    
    if par["var_query_gene_names"]:
        input_modality.var.index = [re.sub("\\.[0-9]+$", "", s) for s in input_modality.var[par["var_query_gene_names"]]]

    if par["model"]:
        if par["refererence"]:
            logger.warning("Both 'model' and 'reference' are provided. 'model' will be used to make cell type predictions, 'reference' will be ignored.")
        model = celltypist.models.Model.load(par["model"])
        logger.info(f"Model features: {model.classifier.features}")
        logger.info(input_modality.var)
    
    elif par["reference"]:
        reference_modality = mu.read_h5mu(par["reference"]).mod[par["modality"]]
        
        # TODO: move to the cross check component
        
        if par["var_reference_gene_names"]:
            reference_modality.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_modality.var[par["var_reference_gene_names"]]]
        
        logger.info("Detecting common vars based on ensembl ids")
        common_ens_ids = list(set(reference_modality.var.index).intersection(set(input_modality.var.index)))
        
        logger.info("  reference n_vars: %i", reference_modality.n_vars)
        logger.info("  input n_vars: %i", input_modality.n_vars)
        logger.info("  intersect n_vars: %i", len(common_ens_ids))
        assert len(common_ens_ids) >= 100, "The intersection of genes is too small."
        
        input_matrix = input_modality.layers[par["reference_layer"]] if par["reference_layer"] else input_modality.X
        reference_matrix = reference_modality.layers[par["reference_layer"]] if par["reference_layer"] else reference_modality.X

        if not check_celltypist_format(input_matrix):
            logger.warning("Input data is not in the reccommended format for CellTypist.")
        if not check_celltypist_format(reference_matrix):
            logger.warning("Reference data is not in the reccommended format for CellTypist.")
        
        labels = reference_modality.obs[par["reference_obs_targets"]]
        model = celltypist.train(reference_matrix,
                                labels=labels,
                                genes=reference_modality.var.index,
                                feature_selection=par["feature_selection"],
                                check_expression=par["check_expression"])
        
    else:
        raise ValueError("Either 'model' or 'reference' has to be provided.")
    
    predictions = celltypist.annotate(input_modality,
                                      model,
                                      majority_voting=par["majority_voting"])
    input_modality.obs[par["output_obs_predictions"]] = predictions.predicted_labels["predicted_labels"]
    input_modality.obs[par["output_obs_probability"]] = predictions.probability_matrix.max(axis=1).values
    
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])
    
if __name__ == '__main__':
    main(par)
