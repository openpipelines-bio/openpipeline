import mudata as mu
import numpy as np
import sys
from pynndescent import PyNNDescentTransformer
from sklearn.neighbors import KNeighborsClassifier

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_obsm_features": None,
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "reference_obsm_features": None,
    "reference_obs_targets": ["cell_type"],
    "output": "foo_distance.h5mu",
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "output_uns_parameters": "labels_transfer",
    "output_compression": None,
    "weights": "distance",
    "n_neighbors": 15
}
meta = {
    "resources_dir": "src/labels_transfer/utils"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from helper import check_arguments, get_reference_features, get_query_features


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


def distances_to_affinities(distances):
    # Apply Gaussian kernel to distances
    stds = np.std(distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)
    distances_tilda = np.exp(-np.true_divide(distances, stds))

    # normalize the distances_tilda
    # if the sum of a row of the distances tilda equals 0,
    # set normalized distances for that row to 1
    # else divide the row values by the value of the sum of the row
    distances_tilda_normalized = np.where(
        np.sum(distances_tilda, axis=1, keepdims=True) == 0,
        1,
        distances_tilda / np.sum(distances_tilda, axis=1, keepdims=True)
    )
    return distances_tilda_normalized


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

if par["input_obsm_distances"] and par["reference_obsm_distances"]:
    logger.info("Using pre-calculated distances for KNN classification as provided in `--input_obsm_distances` and `--reference_obsm_distances`.")
    query_neighbors = q_adata.obsm[par["input_obsm_distances"]]
    reference_neighbors = r_adata.obsm[par["reference_obsm_distances"]]

    if query_neighbors.shape[1] != reference_neighbors.shape[1]:
        raise ValueError("The number of neighbors in the query and reference distance matrices do not match. Make sure both distance matrices contain distances to the reference dataset.")

    # Make sure the number of neighbors present in the distance matrix matches the requested number of neighbors in --n_neighbors
    # Otherwise reduce n_neighbors for KNN
    smallest_neighbor_count = min(
        np.diff(query_neighbors.indptr).min(),
        np.diff(reference_neighbors.indptr).min()
    )
    if smallest_neighbor_count < par["n_neighbors"]:
        logger.warning(f"The number of neighbors in the distance matrices is smaller than the requested number of neighbors in --n_neighbors. Reducing n_neighbors to {smallest_neighbor_count} for KNN Classification")
        par["n_neighbors"] = smallest_neighbor_count

elif par["input_obsm_distances"] or par["reference_obsm_distances"]:
    raise ValueError("Make sure to provide both --input_obsm_distances and --reference_obsm_distances if you want to use a pre-calculated distance matrix for KNN classification.")

elif not par["input_obsm_distances"] and not par["reference_obsm_distances"]:
    logger.info("No pre-calculated distances were provided. Calculating distances using the PyNNDescent algorithm.")
    # Generating training and inference data
    train_X = get_reference_features(r_adata, par, logger)
    inference_X = get_query_features(q_adata, par, logger)

    neighbors_transformer = PyNNDescentTransformer(
        n_neighbors=par["n_neighbors"],
        parallel_batch_queries=True,
    )
    neighbors_transformer.fit(train_X)

    # Square sparse matrix with distances to n neighbors in reference data
    query_neighbors = neighbors_transformer.transform(inference_X)
    reference_neighbors = neighbors_transformer.transform(train_X)

# For each target, train a classifier and predict labels
for obs_tar, obs_pred, obs_proba in zip(par["reference_obs_targets"],  par["output_obs_predictions"], par["output_obs_probability"]):
    logger.info(f"Predicting labels for {obs_tar}")

    weights_dict = {
        "uniform": "uniform",
        "distance": "distance",
        "gaussian": distances_to_affinities
    }

    logger.info(f"Using KNN classifier with {par['weights']} weights")
    train_y = r_adata.obs[obs_tar].to_numpy()
    classifier = KNeighborsClassifier(
        n_neighbors=par["n_neighbors"],
        metric="precomputed",
        weights=weights_dict[par["weights"]]
        )
    classifier.fit(X=reference_neighbors, y=train_y)
    predicted_labels = classifier.predict(query_neighbors)
    probabilities = classifier.predict_proba(query_neighbors).max(axis=1)

    # save_results
    logger.info(f"Saving predictions to {obs_pred} and probabilities to {obs_proba} in obs")
    q_adata.obs[obs_pred] = predicted_labels
    q_adata.obs[obs_proba] = probabilities

logger.info(f"Saving output data to {par['output']}")
q_mdata.mod[par['modality']] = q_adata
q_mdata.write_h5mu(par['output'], compression=par['output_compression'])
