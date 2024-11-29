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
    "reference": None,
    # "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    # "model": None,
    "model": "resources_test/annotation_test_data/svm_model.pkl",
    "reference_obs_target": "cell_ontology_class",
    "reference_var_input": None,
    "input_layer": None,
    "reference_layer": None,
    "max_iter": 5000,
    "c_reg": 1,
    "class_weight": "balanced",
    "output_compression": "gzip",
    "input_var_gene_names": None,
    "reference_var_gene_names": "ensemblid",
    "reference_layer": None,
    "output_obs_prediction": "svm_pred",
    "output_obs_probability": "svm_probability",
    "input_reference_gene_overlap": 100
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from cross_check_genes import cross_check_genes
from subset_vars import subset_vars
from set_var_index import set_var_index

logger = setup_logger()

def main():

    if (not par["model"] and not par["reference"]) or (par["model"] and par["reference"]):
        raise ValueError("Make sure to provide either 'model' or 'reference', but not both.")
    logger.info("Reading input data")
    input_mudata = mu.read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]].copy()
    input_modality = set_var_index(input_modality, par["input_var_gene_names"])

    if par["model"]:
        logger.info("Loading a pre-trained model")
        model = pickle.load(open(par["model"], "rb"))
        if hasattr(model, "_feature_names_in"):
            common_genes = cross_check_genes(input_modality.var.index, model._feature_names_in, par["input_reference_gene_overlap"])
            if not len(common_genes) == len(model._feature_names_in):
                raise ValueError("Input dataset does not contain all model features.")
            input_modality = input_modality[:, common_genes]
            input_matrix = input_modality.layers[par["input_layer"]] if par["input_layer"] else input_modality.X

        else:
            logger.warning("Model does not have feature names saved. Could not check overlap of model's features with query genes.")

    elif par["reference"]:
        logger.info("Reading reference data")

        reference_mudata = mu.read_h5mu(par["reference"])
        reference_modality = reference_mudata.mod[par["modality"]].copy()
        reference_modality = set_var_index(reference_modality, par["reference_var_gene_names"])

        # subset to HVG if required
        if par["reference_var_input"]:
            reference_modality = subset_vars(reference_modality, par["reference_var_input"])

        # Query and input require the exact same features
        common_genes = cross_check_genes(input_modality.var.index, reference_modality.var.index, par["input_reference_gene_overlap"])
        reference_modality = reference_modality[:, common_genes]
        input_modality = input_modality[:, common_genes]

        reference_matrix = reference_modality.layers[par["reference_layer"]] if par["reference_layer"] else reference_modality.X
        input_matrix = input_modality.layers[par["input_layer"]] if par["input_layer"] else input_modality.X

        logger.info("Training a model...")
        labels = reference_modality.obs[par["reference_obs_target"]].to_numpy()
        model = CalibratedClassifierCV(svm.LinearSVC(
            C=par["c_reg"],
            max_iter=par["max_iter"],
            class_weight=par["class_weight"] if not par["class_weight"] == "uniform" else None,
            dual="auto",
        ))
        model.fit(reference_matrix, labels)
        model._feature_names_in = reference_modality.var.index

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