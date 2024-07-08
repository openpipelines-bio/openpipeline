import mudata as mu
import numpy as np
from pynndescent import PyNNDescentTransformer
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_obsm_features": "X_pca",
    "reference": "TS_lung_filtered.h5mu",
    "reference_obsm_features": "X_pca",
    "reference_obs_targets": ["cell_type"],
    "output": "foo.h5mu",
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "output_uns_parameters": "labels_transfer",
    "output_compression": None,
    "weights": "uniform",
    "n_neighbors": 15
}
## VIASH END


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

# Generating training and inference data
## train data
if par["reference_obsm_features"]:
    logger.info(f"Using reference .obsm {par["reference_obsm_features"]} for training")
    train_X = r_adata.obsm[par["reference_obsm_features"]]
else:
    logger.info("Using reference .X as features for training")
    train_X = r_adata.X
    
## inference data
if par["input_obsm_features"]:
    logger.info(f"Using query .obsm {par["input_obsm_features"]} for inference")
    inference_X = q_adata.obsm[par["input_obsm_features"]]
else:
    logger.info("Using query .X as features for inference")
    inference_X = q_adata.X

# Ensure output obs predictions and uncertainties are same length as obs targets
if par["output_obs_predictions"]:
    assert len(par["output_obs_predictions"]) == len(par["reference_obs_targets"]), "output_obs_predictions must be same length as reference_obs_targets"
    output_obs_predictions = par["output_obs_predictions"]
else:
    output_obs_predictions = [x + "_pred" for x in par["reference_obs_targets"]]

if par["output_obs_probability"]:   
    assert len(par["output_obs_probability"]) == len(par["reference_obs_targets"]), "output_obs_probability must be same length as reference_obs_targets"
else:
    output_obs_uncertainties = [x + "_probability" for x in par["reference_obs_targets"]]

# For each target, train a classifier and predict labels
for obs_tar, obs_pred, obs_proba in zip(par["reference_obs_targets"], output_obs_predictions, output_obs_uncertainties):
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
