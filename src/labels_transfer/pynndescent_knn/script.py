import mudata as mu
import numpy as np
from pynndescent import PyNNDescentTransformer
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_subset_scgpt_integrated_leiden_qc.h5mu",
    "modality": "rna",
    "input_obsm_features": "X_PCA"
    "reference": "https://zenodo.org/records/7587774/files/TS_Lung_filtered.h5ad",
    "obsm_input": "",
    "obs_label": "",
    # "output": "output.h5mu",
    # "metric": 'cosine',
    "n_neighbors": 15,
    "knn_weights": "uniform",
    # "modality": "rna",
    # "obsm_input": "X_pca",
    "uns_output": "neighbors",
    "obsp_distances": "distances",
    "obsp_connectivities": "connectivities"
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

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Computing a neighborhood graph on modality %s", mod)
adata = mdata.mod[mod]

# Generate train data and prediction data
logger.info(f"Generate train data and prediction data based on reference index {par["reference_index"]}")
ref_idx = adata.obs["_dataset"] == par["reference_index"]
train_X = adata[ref_idx].obsm[par["obsm_input"]]
train_Y = adata[ref_idx].obs[par["obs_label"]].to_numpy()

# Pipeline instantiation
logger.info(f"Instantiate pipeline of PyNNDescentTransformer and KNeighborsClassifier with n_neighbors {par['n_neighbors']} and knn_weights {par['knn_weights'] and weights par["knn_weights"]}")
knn = make_pipeline(
    PyNNDescentTransformer(
        n_neighbors=par["n_neighbors"],
        parallel_batch_queries=True,
    ),
    KNeighborsClassifier(metric="precomputed", weights=par["knn_weights"]),
)

knn.fit(train_X, train_Y)
knn_pred = knn.predict(adata.obsm[par["layer"]])

# save_results
adata.obs[par["obs_knn_predictions"]] = knn_pred

if par["return_probabilities"]:
    adata.obs[par["obs_knn_probabilities"]] = np.max(
        knn.predict_proba(adata.obsm["X_pca_harmony"]), axis=1
    )
