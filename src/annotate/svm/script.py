import sys
import logging
import mudata as mu
import anndata as ad
import re
import numpy as np
import os
from tqdm import tqdm
from sklearn.calibration import CalibratedClassifierCV
from sklearn import svm
import pickle



## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "model": None,
    "reference_obs_targets": ["cell_ontology_class"],
    "input_layer": None,
    "reference_layer": None,
    "max_iter": 1000,
    "c_reg": 1,
    "class_weight": "balanced",
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "output_compression": "gzip"
}
meta = {"resources_dir": "src/annotate/svm"}
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

logger = setup_logger()

def main():
    logger.info("Reading input data")
    input_mudata = mu.read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]].copy()

    obs_predictions = par["output_obs_predictions"] if par["output_obs_predictions"] else [f"{target}_pred" for target in par["reference_obs_targets"]]
    obs_probabilities = par["output_obs_probability"] if par["output_obs_probability"] else [f"{target}_prob" for target in par["reference_obs_targets"]]

    input_matrix = input_modality.layers[par["input_layer"]] if par["input_layer"] else input_modality.X 

    models = {}

    if par["reference"]:
        logger.info("Reading reference data")

        reference_mudata = mu.read_h5mu(par["reference"])
        reference_modality = reference_mudata.mod[par["modality"]].copy()

        reference_modality.var["gene_symbol"] = list(reference_modality.var.index)
        reference_modality.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_modality.var["ensemblid"]]

        logger.info("Detecting common vars based on ensembl ids")
        common_ens_ids = list(set(reference_modality.var.index).intersection(set(input_modality.var.index)))

        logger.info("  reference n_vars: %i", reference_modality.n_vars)
        logger.info("  input n_vars: %i", input_modality.n_vars)
        logger.info("  intersect n_vars: %i", len(common_ens_ids))
        assert len(common_ens_ids) >= 100, "The intersection of genes is too small."

        reference_matrix = reference_modality.layers[par["reference_layer"]] if par["reference_layer"] else reference_modality.X

        logger.info("Training models...")
        for reference_obs_target in tqdm(par["reference_obs_targets"]):

            logger.info(f"Training model for {reference_obs_target}")
            labels = reference_modality.obs[reference_obs_target].to_numpy()
            model = CalibratedClassifierCV(svm.LinearSVC(
                C=par["c_reg"],
                max_iter=par["max_iter"],
                class_weight=par["class_weight"],
                verbose=1
            ))
            model.fit(reference_matrix, labels)
            models[reference_obs_target] = model

    elif par["model"]:
        logger.info("Loading pre-trained models")
        for model_path, reference_obs_target  in zip(par["model"], par["reference_obs_targets"]):
            logger.info(f"Loading model for {reference_obs_target}")
            model = pickle.load(open(model_path, "rb"))
            models[reference_obs_target] = model

    else:
        raise ValueError("Either reference or model must be provided")
    
    logger.info("Running predictions...")
    for model, obs_prediction, obs_probability in models.values(), obs_predictions, obs_probabilities:
        predictions = model.predict(input_matrix)
        probabilities = np.max(model.predict_proba(input_modality), axis=1)
        
        input_modality.obs[obs_prediction] = predictions
        input_modality.obs[obs_probability] = probabilities

    logger.info("Writing output data")
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])

if __name__ == "__main__":
    main()