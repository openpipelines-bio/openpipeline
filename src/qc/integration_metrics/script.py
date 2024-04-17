import scib
import scanpy as sc
import mudata as mu
import numpy as np


## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_integrated.h5mu",
    "modality": "rna",
    "var_gene_names": None,
    "uns_neighbors": "scGPT_integration_neighbors",
    "obsp_neighbor_distances": "scGPT_integration_distances",
    "obsp_neighbor_connectivities": "scGPT_integration_connectivities",
    "obsm_embeddings": "X_scGPT",
    "obsm_umap": "X_scGPT_umap",
    "obs_batch_label": "sample",
    "obs_cell_label": "cell_type",
    "output_compression": None,
    "obs_cluster": "louvain_cluster",
    "seed": None,
    "output": "resources_test/scgpt/test_resources/Kim2020_Lung_integrated_qc.h5mu",
}

## VIASH END
if par["seed"]:
    np.random.seed(par["seed"])

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
logger = setup_logger()

logger.info("Reading in data")
# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

# Ensure all batch labels are str
adata.obs[par["obs_batch_label"]] = adata.obs[par["obs_batch_label"]].astype(str)

# Remove nan from cell labels
if adata.obs["cell_type"].isna().any():
    logger.warning("Cell label obs column contains nan values: removing rows with nan value")
    adata_t = adata[~adata.obs["cell_type"].isna()].copy()
else:
    adata_t = adata.copy()

if par["obs_cell_label"]:
    logger.info("Generating clusters")
    res_max, nmi_max, nmi_all = scib.metrics.clustering.cluster_optimal_resolution(
                adata_t,
                label_key=par["obs_cell_label"],
                cluster_key=par["obs_cluster"],
                cluster_function=sc.tl.louvain,
                use_rep=par["obsm_embeddings"],
                metric=scib.metrics.nmi,
                verbose=False,
                return_all=True,
                force=True,
            )

    logger.info("Calculating NMI score")
    nmi_score = scib.metrics.nmi(
            adata_t,
            cluster_key=par["obs_cluster"],
            label_key=par["obs_cell_label"],
            implementation="arithmetic",
            nmi_dir=None,
        )

    logger.info("Calculating ARI score")
    ari_score = scib.metrics.ari(
        adata_t,
        cluster_key=par["obs_cluster"],
        label_key=par["obs_cell_label"]
    )

    logger.info("Calculating AWS score (cell types)")
    asw_label = scib.metrics.silhouette(
        adata_t,
        label_key=par["obs_cell_label"],
        embed=par["obsm_embeddings"],
        metric="euclidean"
        )

    logger.info("Calculating AWS score (batches)")
    asw_batch = scib.metrics.silhouette_batch(
            adata_t,
            batch_key=par["obs_batch_label"],
            label_key=par["obs_cell_label"],
            embed=par["obsm_embeddings"],
            metric="euclidean",
            return_all=False,
            verbose=False,
        )

    logger.info("Calculating PC regression")
    pcr_score = scib.metrics.pcr_comparison(
        adata_t,
        adata_t,
        embed=par["obsm_embeddings"],
        covariate=par["obs_batch_label"],
        verbose=False
    )

    logger.info("Calculating graph connectivity")
    graph_conn_score = scib.metrics.graph_connectivity(
        adata_t,
        label_key=par["obs_cell_label"]
    )

logger.info("Calculating average bio score")

avg_bio = np.mean(
    [
        nmi_score, ari_score, asw_label
    ]
)

logger.info("Writing output data")
adata.uns["ari_score"] = ari_score
adata.uns["nmi_score"] = nmi_score
adata.uns["asw_label"] = asw_label
adata.uns["asw_batch"] = asw_batch
adata.uns["pcr_score"] = pcr_score
adata.uns["graph_conn_score"] = graph_conn_score
adata.uns["avg_bio"] = avg_bio

mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])