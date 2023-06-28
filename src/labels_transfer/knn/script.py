import logging
import warnings

import mudata
import numpy as np
import scanpy as sc
from scipy.sparse import issparse
import pynndescent
import numba


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "reference": "https://zenodo.org/record/6337966/files/HLCA_emb_and_metadata.h5ad",
    "targets": ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"],
    "modality": "rna",
    "reference_obsm_key": "X_integrated_scanvi",
    "query_obsm_key": "X_integrated_scanvi",
    "output": "foo.h5mu",
    "output_obs_suffix": "_pred",
    "output_uns_key": "labels_transfer",
    "n_neighbors": 1
}
### VIASH END

def _setup_logger():
    # TODO: move to utils 
    
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger


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

def main():
    logger = _setup_logger()

    logger.info("Reading input (query) data")
    mdata = mudata.read(par["input"].strip())
    adata = mdata.mod[par["modality"]]

    logger.info("Reading reference data")
    adata_reference = sc.read(par["reference"], backup_url=par["reference"])

    if par["reference_obsm_key"] is None:
        logger.info("Using X layer of reference data")
        X_train = adata_reference.X
    else:
        logger.info(f"Using obsm layer {par['reference_obsm_key']} of reference data")
        X_train = adata_reference.obsm[par["reference_obsm_key"]]

    # pynndescent does not support sparse matrices
    if issparse(X_train):
        warnings.warn("Converting sparse matrix to dense. This may consume a lot of memory.")
        X_train = X_train.toarray()

    logger.debug(f"Shape of train data: {X_train.shape}")

    logger.info("Building NN index")
    ref_nn_index = pynndescent.NNDescent(X_train, n_neighbors=par["n_neighbors"])
    ref_nn_index.prepare()

    if par["query_obsm_key"] is None:
        logger.info("Using X layer of query data")
        query = adata.X
    else:
        logger.info(f"Using obsm layer {par['query_obsm_key']} of query data")
        query = adata.obsm[par["query_obsm_key"]]

    ref_neighbors, ref_distances = ref_nn_index.query(query, k=par["n_neighbors"])

    weights = distances_to_affinities(ref_distances)

    if par["output_uns_key"] not in adata.uns:
        adata.uns[par["output_uns_key"]] = {}

    # for each annotation level, get prediction and uncertainty
    for target in par["targets"]:
        logger.info(f"Predicting labels for {target}")
        ref_cats = adata_reference.obs[target].cat.codes.to_numpy()[ref_neighbors]
        prediction, uncertainty = weighted_prediction(weights, ref_cats)
        prediction = np.asarray(adata_reference.obs[target].cat.categories)[prediction]
        
        predicted_label_col_name = target + par["output_obs_suffix"]
        adata.obs[predicted_label_col_name], adata.obs[target + "_uncertainty"] = prediction, uncertainty
        
        # Write information about labels transfer to uns
        adata.uns[par["output_uns_key"]][predicted_label_col_name] = {
            "method": "KNN_pynndescent",
            "n_neighbors": par["n_neighbors"],
            "reference": par["reference"]
        }

    mdata.mod[par['modality']] = adata
    mdata.update()
    mdata.write_h5mu(par['output'].strip())

if __name__ == "__main__":
    main()
