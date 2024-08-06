import mudata as mu
import numba
import numpy as np
from scipy.sparse import issparse
import sys
from pynndescent import PyNNDescentTransformer
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_obsm_features": None,
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "reference_obsm_features": None,
    "reference_obs_targets": ["cell_type"],
    "output": "foo.h5mu",
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "output_uns_parameters": "labels_transfer",
    "output_compression": None,
    "weights": "uniform",
    "n_neighbors": 15
}
meta ={
    "resources_dir": "src/labels_transfer/utils"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from helper import check_arguments, get_reference_features, get_query_features, check_sparsity

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

@numba.njit
def weighted_prediction(weights, ref_cats):
    """Get highest weight category."""
    N = len(weights)
    predictions = np.zeros((N,), dtype=ref_cats.dtype)
    probabilities = np.zeros((N,))
    for i in range(N):
        obs_weights = weights[i]
        obs_cats = ref_cats[i]
        best_prob = 0
        for c in np.unique(obs_cats):
            cand_prob = np.sum(obs_weights[obs_cats == c])
            if cand_prob > best_prob:
                best_prob = cand_prob
                predictions[i] = c
                probabilities[i] = best_prob

    return predictions, probabilities

def distances_to_affinities(distances):
    stds = np.std(distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)
    distances_tilda = np.exp(-np.true_divide(distances, stds))

    return distances_tilda / np.sum(distances_tilda, axis=1, keepdims=True)

logger = setup_logger()

# Reading in data
logger.info(f"Reading in query dataset {par['input']} and reference datasets {par['reference']}")
q_mdata = mu.read_h5mu(par["input"])
q_adata = q_mdata.mod[par["modality"]]

r_mdata = mu.read_h5mu(par["reference"])
r_adata = r_mdata.mod[par["modality"]]

# check arguments
logger.info("Checking arguments")
par = check_arguments(par)

# Generating training and inference data
logger.info("Generating training and inference data")
train_X = get_reference_features(r_adata, par, logger)
inference_X = get_query_features(q_adata, par, logger)

# pynndescent does not support sparse matrices
train_X = check_sparsity(train_X, logger)
inference_X = check_sparsity(inference_X, logger)

neighbors_transformer = PyNNDescentTransformer(
    n_neighbors=par["n_neighbors"],
    parallel_batch_queries=True,
)
neighbors_transformer.fit(train_X)

# Square sparse matrix with distances to n neighbors in reference data
reference_neighbors = neighbors_transformer.transform(inference_X)

# For each target, train a classifier and predict labels
for obs_tar, obs_pred, obs_proba in zip(par["reference_obs_targets"],  par["output_obs_predictions"], par["output_obs_probability"]):
    logger.info(f"Predicting labels for {obs_tar}")

    if par["weights"] != "gaussian":
        train_y = r_adata.obs[obs_tar].to_numpy()
        classifier = KNeighborsClassifier(n_neighbors=50, metric="precomputed", weights=par["weights"])
        classifier.fit(
            X=neighbors_transformer.transform(train_X), y=train_y
        )
        predicted_labels = classifier.predict(reference_neighbors)
        probabilities = classifier.predict_proba(reference_neighbors).max(axis=1)
        
    elif par["weights"] == "gaussian":
        # Convert type to category so that the code below works for any initial type
        r_adata.obs[obs_tar] = r_adata.obs[obs_tar].astype("category")

        ref_distances = reference_neighbors.data.reshape(inference_X.shape[0], par["n_neighbors"])
        ref_neighbors_idxs = reference_neighbors.indices.reshape(inference_X.shape[0], par["n_neighbors"])

        # Get indices of each category because numba throws an error when using strings
        ref_cats = r_adata.obs[obs_tar].cat.codes.to_numpy()[ref_neighbors_idxs]
        affinities = distances_to_affinities(ref_distances)
        predicted_labels, probabilities = weighted_prediction(affinities, ref_cats)
        # Convert indices back to readable labels
        predicted_labels = np.asarray(r_adata.obs[obs_tar].cat.categories)[predicted_labels]
        
    else:
        raise ValueError("Wriong weights parameter")

    logger.info(f"Predicting {obs_pred} predictions and {obs_proba} probabilities")

    # save_results
    q_adata.obs[obs_pred] = predicted_labels
    q_adata.obs[obs_proba] = probabilities

logger.info(f"Saving output data to {par['output']}")
q_mdata.mod[par['modality']] = q_adata
q_mdata.write_h5mu(par['output'], compression=par['output_compression'])
