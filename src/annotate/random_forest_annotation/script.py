import sys
import mudata as mu
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import pickle


## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "input_var_gene_names": None,
    "input_reference_gene_overlap": 100,
    "modality": "rna",
    "reference": None,
    # "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    # "model": None,
    "model": "resources_test/annotation_test_data/rf_model.pkl",
    "reference_obs_target": "cell_ontology_class",
    "reference_var_gene_names": "ensemblid",
    "reference_var_input": None,
    "input_layer": None,
    "reference_layer": None,
    "n_estimators": 100,
    "criterion": "gini",
    "max_depth": None,
    "class_weight": None,
    "max_features": 200,
    "output_compression": "gzip",
    "output_obs_predictions": "random_forest_pred",
    "output_obs_probability": "random_forest_probability",
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
    logger.info("Reading input data")
    input_mudata = mu.read_h5mu(par["input"])
    input_adata = input_mudata.mod[par["modality"]]
    input_modality = input_adata.copy()
    input_modality = set_var_index(
        input_modality, par["input_var_gene_names"], par["sanitize_gene_names"]
    )

    # Handle max_features parameter
    max_features_conversion = {
        "all": None,
        "sqrt": "sqrt",
        "log2": "log2",
    }
    try:
        max_features = max_features_conversion.get(
            par["max_features"], int(par["max_features"])
        )
    except ValueError:
        raise ValueError(
            f"Invaldid value {par['max_features']} for --max_features: must either be an integer or one of 'sqrt', 'log2' or 'all'"
        )

    if (not par["model"] and not par["reference"]) or (
        par["model"] and par["reference"]
    ):
        raise ValueError(
            "Make sure to provide either 'model' or 'reference', but not both."
        )

    if par["model"]:
        logger.info("Loading a pre-trained model")
        model = pickle.load(open(par["model"], "rb"))
        if hasattr(model, "_feature_names_in"):
            common_genes = cross_check_genes(
                input_modality.var.index,
                model._feature_names_in,
                par["input_reference_gene_overlap"],
            )
            if not len(common_genes) == len(model._feature_names_in):
                raise ValueError("Input dataset does not contain all model features.")
            input_modality = input_modality[:, common_genes]
            input_matrix = (
                input_modality.layers[par["input_layer"]]
                if par["input_layer"]
                else input_modality.X
            )

        else:
            logger.warning(
                "Model does not have feature names saved. Could not check overlap of model's features with query genes."
            )

    elif par["reference"]:
        logger.info("Reading reference data")

        reference_mudata = mu.read_h5mu(par["reference"])
        reference_modality = reference_mudata.mod[par["modality"]].copy()
        reference_modality = set_var_index(
            reference_modality,
            par["reference_var_gene_names"],
            par["sanitize_gene_names"],
        )

        # subset to HVG if required
        if par["reference_var_input"]:
            reference_modality = subset_vars(
                reference_modality, par["reference_var_input"]
            )

        # Query and input require the exact same features
        common_genes = cross_check_genes(
            input_modality.var.index,
            reference_modality.var.index,
            par["input_reference_gene_overlap"],
        )
        reference_modality = reference_modality[:, common_genes]
        input_modality = input_modality[:, common_genes]

        reference_matrix = (
            reference_modality.layers[par["reference_layer"]]
            if par["reference_layer"]
            else reference_modality.X
        )
        input_matrix = (
            input_modality.layers[par["input_layer"]]
            if par["input_layer"]
            else input_modality.X
        )

        logger.info("Training a model...")
        labels = reference_modality.obs[par["reference_obs_target"]].to_numpy()
        model = RandomForestClassifier(
            n_estimators=par["n_estimators"],
            criterion=par["criterion"],
            max_depth=par["max_depth"],
            class_weight=par["class_weight"]
            if not par["class_weight"] == "uniform"
            else None,
            max_features=max_features,
        )
        model.fit(reference_matrix, labels)
        model._feature_names_in = reference_modality.var.index

    logger.info("Running predictions...")
    predictions = model.predict(input_matrix)
    probabilities = np.max(model.predict_proba(input_matrix), axis=1)

    input_adata.obs[par["output_obs_predictions"]] = predictions
    input_adata.obs[par["output_obs_probability"]] = probabilities

    logger.info("Writing output data")
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
