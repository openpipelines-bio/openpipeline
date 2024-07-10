import mudata as mu
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
from helper_2 import check_arguments, get_reference_features, get_query_features, check_sparsity

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
        
# For each target, train a classifier and predict labels
for obs_tar, obs_pred, obs_proba in zip(par["reference_obs_targets"],  par["output_obs_predictions"], par["output_obs_probability"]):
    logger.info(f"Predicting labels for {obs_tar}")
    train_Y = r_adata.obs[obs_tar].to_numpy()

    # Pipeline instantiation
    logger.info(f"Instantiate pipeline of PyNNDescentTransformer and KNeighborClassifier with {par['n_neighbors']} n_neighbors and {par['weights']} weights")
    knn = make_pipeline(
        PyNNDescentTransformer(
            n_neighbors=par["n_neighbors"],
            parallel_batch_queries=True,
        ),
        KNeighborsClassifier(metric="precomputed", weights=par["weights"]),
    )

    logger.info(f"Training PyNNDescentTransformer and KNeighborClassifier based on {obs_tar} obs labels")
    knn.fit(train_X, train_Y)

    logger.info(f"Predicting {obs_pred} predictions and {obs_proba} probabilities")
    knn_pred = knn.predict(inference_X)
    knn_proba = knn.predict_proba(inference_X)

    # save_results
    q_adata.obs[obs_pred] = knn_pred
    q_adata.obs[obs_proba] = np.max(knn_proba, axis=1)

logger.info(f"Saving output data to {par['output']}")
q_mdata.mod[par['modality']] = q_adata
q_mdata.write_h5mu(par['output'], compression=par['output_compression'])
