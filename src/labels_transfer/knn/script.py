import sys
import warnings

import mudata
import numpy as np
import scanpy as sc
from scipy.sparse import issparse
import pynndescent
import numba


## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_obsm_features": "X_pca",
    "reference": "https://zenodo.org/records/7587774/files/TS_Blood_filtered.h5ad",
    "reference_obsm_features": "X_pca",
    "reference_obs_targets": ["cell_type"],
    "output": "foo.h5mu",
    "output_obs_predictions": None,
    "output_obs_uncertainty": None,
    "output_uns_parameters": "labels_transfer",
    "n_neighbors": 1
}
meta = {
    "resources_dir": "src/labels_transfer/utils"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from helper import check_arguments, get_reference_features, get_query_features
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
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
    uncertainty = np.zeros((N,))
    for i in range(N):
        obs_weights = weights[i]
        obs_cats = ref_cats[i]
        best_prob = 0
        for c in np.unique(obs_cats):
            cand_prob = np.sum(obs_weights[obs_cats == c])
            if cand_prob > best_prob:
                best_prob = cand_prob
                predictions[i] = c
                uncertainty[i] = max(1 - best_prob, 0)

    return predictions, uncertainty

def distances_to_affinities(distances):
    stds = np.std(distances, axis=1)
    stds = (2.0 / stds) ** 2
    stds = stds.reshape(-1, 1)
    distances_tilda = np.exp(-np.true_divide(distances, stds))

    return distances_tilda / np.sum(distances_tilda, axis=1, keepdims=True)

def main(par):
    logger = setup_logger()

    logger.info("Checking arguments")
    par = check_arguments(par)

    logger.info("Reading input (query) data")
    mdata = mudata.read(par["input"])
    adata = mdata.mod[par["modality"]]

    logger.info("Reading reference data")
    adata_reference = sc.read(par["reference"], backup_url=par["reference"])

    # fetch feature data
    train_data = get_reference_features(adata_reference, par, logger)
    query_data = get_query_features(adata, par, logger)

    # pynndescent does not support sparse matrices
    if issparse(train_data):
        warnings.warn("Converting sparse matrix to dense. This may consume a lot of memory.")
        train_data = train_data.toarray()

    logger.debug(f"Shape of train data: {train_data.shape}")

    logger.info("Building NN index")
    ref_nn_index = pynndescent.NNDescent(train_data, n_neighbors=par["n_neighbors"])
    ref_nn_index.prepare()

    ref_neighbors, ref_distances = ref_nn_index.query(query_data, k=par["n_neighbors"])

    weights = distances_to_affinities(ref_distances)

    output_uns_parameters = adata.uns.get(par["output_uns_parameters"], {})

    # for each annotation level, get prediction and uncertainty
    
    for obs_tar, obs_pred, obs_unc in zip(par["reference_obs_targets"], par["output_obs_predictions"], par["output_obs_uncertainty"]):
        logger.info(f"Predicting labels for {obs_tar}")
        ref_cats = adata_reference.obs[obs_tar].cat.codes.to_numpy()[ref_neighbors]
        prediction, uncertainty = weighted_prediction(weights, ref_cats)
        prediction = np.asarray(adata_reference.obs[obs_tar].cat.categories)[prediction]
        
        adata.obs[obs_pred], adata.obs[obs_unc] = prediction, uncertainty
        
        # Write information about labels transfer to uns
        output_uns_parameters[obs_tar] = {
            "method": "KNN_pynndescent",
            "n_neighbors": par["n_neighbors"],
            "reference": par["reference"]
        }

    adata.uns[par["output_uns_parameters"]] = output_uns_parameters

    mdata.mod[par['modality']] = adata
    mdata.update()
    mdata.write_h5mu(par['output'].strip())

if __name__ == "__main__":
    main(par)
