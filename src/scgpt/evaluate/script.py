from typing import List
import anndata as ad
import scanpy as sc
import scib
from scgpt.utils import eval_scib_metrics

## VIASH START
par = {
    "input": "src/scgpt/test_resources/Kim2020_Lung_embedded.h5ad",
    "embedding_layer": "X_scGPT",
    "cell_type_layer": "cell_type",
    "batch_layer": "sample",
    "umap_embeddings_batch": "Kim2020_Lung_batch_umap.png",
    "umap_embeddings_celltype": "Kim2020_Lung_celltype_umap.png",
}
## VIASH END

input_adata = ad.read_h5ad(par["input"])
adata = input_adata.copy()

# # Generate umaps
# sc.pp.neighbors(adata, use_rep=par["embedding_layer"])
# sc.tl.umap(adata)
# batch_umap_fig = sc.pl.umap(
#     adata,
#     color=[par["batch_layer"]],
#     frameon=False,
#     return_fig=True,
#     show=False,
#     title="scGPT zero-shot: batch label"
# )
# celltype_umap_fig = sc.pl.umap(
#     adata,
#     color=[par["cell_type_layer"]],
#     frameon=False,
#     return_fig=True,
#     show=False,
#     title="scGPT zero-shot: cell type"
# )

# # Save umap figs
# batch_umap_fig.savefig(
#     par["umap_embeddings_batch"], dpi=300
# )

# celltype_umap_fig.savefig(
#     par["umap_embeddings_celltype"], dpi=300
# )

# Calculate eval metrics
eval_results = scib.metrics.metrics(
    adata,
    adata_int=adata,
    batch_key=par["batch_layer"],
    label_key=par["cell_type_layer"],
    embed=par["embedding_layer"],
    silhouette_=True,
    hvg_score_=False,
    graph_conn_=True,
    pcr_=True,
    isolated_labels_f1_=False,
    trajectory_=False,
    nmi_=True,  # use the clustering, bias to the best matching
    ari_=True,  # use the clustering, bias to the best matching
    cell_cycle_=False,
    ilisi_=False,
    clisi_=False
    )

print(eval_results)
