import sys
import logging
import mudata as mu
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import pickle


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
    "n_estimators": 100,
    "criterion": "gini",
    "max_depth": None,
    "class_weight": None,
    "max_features": 200,
    "output_compression": "gzip",
    "reference_layer": None,
    "output_obs_predictions": "random_forest_pred",
    "output_obs_probability": "random_forest_probability"
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

    input_matrix = input_modality.layers[par["input_layer"]] if par["input_layer"] else input_modality.X 
    
    # Handle max_features parameter
    max_features_conversion = {
        "all": None,
        "sqrt": "sqrt",
        "log2": "log2",  
    }
    try:
        max_features = max_features_conversion.get(par["max_features"], int(par["max_features"]))
    except ValueError:
        raise ValueError(f"Invaldid value {par['max_features']} for --max_features: must either be an integer or one of \'sqrt\', \'log2\' or \'all\'")
        
    if (not par["model"] and not par["reference"]) or (par["model"] and par["reference"]):
        raise ValueError("Make sure to provide either 'model' or 'reference', but not both.")
    
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
        model = RandomForestClassifier(
            n_estimators=par["n_estimators"],
            criterion=par["criterion"],
            max_depth=par["max_depth"],
            class_weight=par["class_weight"] if not par["class_weight"] == "uniform" else None,
            max_features=max_features
        )
        model.fit(reference_matrix, labels)

    logger.info("Running predictions...")
    predictions = model.predict(input_matrix)
    probabilities = np.max(model.predict_proba(input_matrix), axis=1)

    input_modality.obs[par["output_obs_predictions"]] = predictions
    input_modality.obs[par["output_obs_probability"]] = probabilities

    logger.info("Writing output data")
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])

if __name__ == "__main__":
    main()