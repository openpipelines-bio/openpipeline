import sys
import logging
import celltypist
import mudata as mu
import anndata as ad
import re
import numpy as np

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/tmp_TS_Blood_filtered.h5ad",
    "model": None,
    "reference_obs_targets": "cell_ontology_class",
    "check_expression": False,
    "feature_selection": True,
    "majority_voting": True,
    "output_compression": "gzip"
}
meta = {}
## VIASH END

sys.path.append(meta["resources_dir"])
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
        
    if par["model"]:
        model = celltypist.models.Model.load(par["model"])
    
    elif par["reference"]:
        reference_adata = mu.read_h5mu(par["reference"]).mod[par["modality"]]
        
        reference_adata.var["gene_symbol"] = list(reference_adata.var.index)
        reference_adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_adata.var["ensemblid"]]
        
        logger.info("Detecting common vars based on ensembl ids")
        common_ens_ids = list(set(reference_adata.var.index).intersection(set(input_modality.var.index)))
        
        logger.info("  reference n_vars: %i", reference_adata.n_vars)
        logger.info("  input n_vars: %i", input_modality.n_vars)
        logger.info("  intersect n_vars: %i", len(common_ens_ids))
        assert len(common_ens_ids) >= 100, "The intersection of genes is too small."
        
        if par["input_layer"]:
            input_matrix = input_modality.layers[par["input_layer"]]
        else:
            input_matrix = input_modality.X
        if par["reference_layer"]:
            reference_matrix = reference_adata.layers[par["reference_layer"]]
        else:
            reference_matrix = reference_adata.X
        
        if not check_celltypist_format(input_matrix):
            logger.warning("Input data is not in the reccommended format for CellTypist.")
        if not check_celltypist_format(reference_matrix):
            logger.warning("Reference data is not in the reccommended format for CellTypist.")
        
        models = []
        for label in par["reference_obs_targets"]:
            labels = reference_adata.obs[label]
            model = celltypist.train(reference_matrix,
                                    labels = labels,
                                    genes = reference_adata.var_names,
                                    feature_selection=par["feature_selection"],
                                    check_expression=par["check_expression"])
            models.append(model)
        
    else:
        raise ValueError("Either 'model' or 'reference' has to be provided.")
    
    obs_predictions = par["output_obs_predictions"] if par["output_obs_predictions"] else ["pred"] * len(models)
    obs_probabilities = par["output_obs_probability"] if par["output_obs_probability"] else ["prob"] * len(models)
    
    for model, reference_obs_target, obs_prediction, obs_probability in zip(models, par["reference_obs_targets"], obs_predictions, obs_probabilities):
        predictions = celltypist.annotate(input_modality,
                                          model,
                                          majority_voting=par["majority_voting"])
        input_modality.obs[f"{reference_obs_target}_{obs_prediction}"] = predictions.predicted_labels["predicted_labels"]
        input_modality.obs[f"{reference_obs_target}_{obs_probability}"] = predictions.probability_matrix.max(axis=1).values
    
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])
    
if __name__ == '__main__':
    main(par)
