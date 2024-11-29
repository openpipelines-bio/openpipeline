import sys
import logging
import mudata as mu
import numpy as np
from sklearn.calibration import CalibratedClassifierCV
from sklearn import svm
import pickle
import re


## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "model": None,
    "reference_obs_target": "cell_ontology_class",
    "input_layer": None,
    "reference_layer": None,
    "max_iter": 5000,
    "c_reg": 1,
    "class_weight": "balanced",
    "output_compression": "gzip",
    "var_query_gene_names": None,
    "var_reference_gene_names": "ensemblid",
    "reference_layer": None,
    "output_obs_prediction": "svm_pred",
    "output_obs_probability": "svm_probability",
}
meta = {"resources_dir": "src/annotate/svm"}
## VIASH END

sys.path.append(meta["resources_dir"])


logger = setup_logger()

def main():
    
    if (not par["model"] and not par["reference"]) or (par["model"] and par["reference"]):
        raise ValueError("Make sure to provide either 'model' or 'reference', but not both.")
    
    logger.info("Reading input data")
    input_mudata = mu.read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]].copy()
    
    input_matrix = input_modality.layers[par["input_layer"]] if par["input_layer"] else input_modality.X 
    
    if par["model"]:
        logger.info("Loading a pre-trained model")
        model = pickle.load(open(par["model"], "rb"))

    elif par["reference"]:
        logger.info("Reading reference data")
        
        reference_mudata = mu.read_h5mu(par["reference"])
        reference_modality = reference_mudata.mod[par["modality"]].copy()
    
        reference_matrix = reference_modality.layers[par["reference_layer"]] if par["reference_layer"] else reference_modality.X

        logger.info("Training a model...")
        labels = reference_modality.obs[par["reference_obs_target"]].to_numpy()
        model = CalibratedClassifierCV(svm.LinearSVC(
            C=par["c_reg"],
            max_iter=par["max_iter"],
            class_weight=par["class_weight"] if not par["class_weight"]=="uniform" else None,
            dual="auto",
        ))
        model.fit(reference_matrix, labels)
    
    logger.info("Running predictions...")
    predictions = model.predict(input_matrix)
    probabilities = np.max(model.predict_proba(input_matrix), axis=1)
    
    input_modality.obs[par["output_obs_prediction"]] = predictions
    input_modality.obs[par["output_obs_probability"]] = probabilities

    logger.info("Writing output data")
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])

if __name__ == "__main__":
    main()