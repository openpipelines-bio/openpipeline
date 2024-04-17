import scib
import scanpy as sc
import mudata as mu
import numpy as np


## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_integrated_qc.h5mu",
    "modality": "rna",
    # "input_layer": None,
    "var_gene_names": None,
    "uns_neighbors": "scGPT_integration_neighbors",
    # "obsp_neighbor_distances": "scGPT_integration_distances",
    # "obsp_neighbor_connectivities": "scGPT_integration_connectivities",
    # "obsm_embeddings": "X_scGPT",
    "obsm_umap": "X_scGPT_umap",
    "obs_batch_label": "sample",
    "obs_cell_label": "cell_type",
    # "output_compression": None,
    # "obs_cluster": "louvain_cluster",
    # "seed": 0,
    "output": "output.qmd",
}

## VIASH END

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

# Open template
with open('src/qc/integration_qc_report/template.md', 'r') as file:
    template = file.read()

logger.info("Generating figures")
logger.info("Generating plots")
fig_cell_type_clusters = sc.pl.embedding(
    adata,
    basis=par["obsm_umap"],
    gene_symbols=par["var_gene_names"],
    neighbors_key=par["uns_neighbors"],
    color=par["obs_cell_label"],
    title="UMAP visualization (label)",
    frameon=False,
    return_fig=True,
    show=False,
)
# fig_cell_type_clusters.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

fig_batch_clusters = sc.pl.embedding(
    adata,
    basis=par["obsm_umap"],
    gene_symbols=par["var_gene_names"],
    neighbors_key=par["uns_neighbors"],
    color=par["obs_batch_label"],
    title="UMAP visualization (batch)",
    frameon=False,
    return_fig=True,
    show=False,
)
# fig_batch_clusters.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)



fig_cell_type_clusters.savefig("umap_label.png", dpi=300, bbox_inches='tight')
fig_batch_clusters.savefig("umap_batch.png", dpi=300, bbox_inches='tight')

# Generate variables required for report
variables = {
    "ari_score": round(adata.uns["ari_score"], 4),
    "nmi_score": round(adata.uns["nmi_score"], 4),
    "asw_label": round(adata.uns["asw_label"], 4),
    "asw_batch": round(adata.uns["asw_batch"], 4),
    "pcr_score": round(adata.uns["pcr_score"], 4),
    "graph_conn_score": round(adata.uns["graph_conn_score"], 4),
    "avg_bio": round(adata.uns["avg_bio"], 4),
    "umap_batch": "umap_batch.png",
    "umap_label": "umap_label.png",
}

markdown_content = template.format(**variables)

with open(par["output"], 'w') as f:
    f.write(markdown_content)
