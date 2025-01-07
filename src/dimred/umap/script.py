import scanpy as sc
import mudata as mu
import sys
import anndata as ad

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "obsm_output": "X_umap",
    "min_dist": 0.5,
    "spread": 1.0,
    "num_components": 2,
    "max_iter": None,
    "alpha": 1.0,
    "gamma": 1.0,
    "negative_sample_rate": 5,
    "init_pos": "spectral",
    "uns_neighbors": "neighbors",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading %s, modality %s", par["input"], par["modality"])
data = mu.read_h5ad(par["input"], mod=par["modality"])

logger.info("Computing UMAP for modality '%s'", par["modality"])

if par["uns_neighbors"] not in data.uns:
    raise ValueError(
        f"'{par['uns_neighbors']}' was not found in .mod['{par['modality']}'].uns."
    )

# create temporary AnnData
# ... because sc.tl.umap doesn't allow to choose
# the obsm output slot
# ... also we can see scanpy is a data format dependency hell
neigh_key = par["uns_neighbors"]
temp_uns = {neigh_key: data.uns[neigh_key]}
conn_key = temp_uns[neigh_key]["connectivities_key"]
dist_key = temp_uns[neigh_key]["distances_key"]
temp_obsp = {
    conn_key: data.obsp[conn_key],
    dist_key: data.obsp[dist_key],
}
pca_key = temp_uns[neigh_key]["params"]["use_rep"]
temp_obsm = {pca_key: data.obsm[pca_key]}

temp_adata = ad.AnnData(obsm=temp_obsm, obsp=temp_obsp, uns=temp_uns, shape=data.shape)

sc.tl.umap(
    temp_adata,
    min_dist=par["min_dist"],
    spread=par["spread"],
    n_components=par["num_components"],
    maxiter=par["max_iter"],
    alpha=par["alpha"],
    gamma=par["gamma"],
    negative_sample_rate=par["negative_sample_rate"],
    init_pos=par["init_pos"],
    neighbors_key=neigh_key,
)

data.obsm[par["obsm_output"]] = temp_adata.obsm["X_umap"]

logger.info("Writing to %s, compression %s", par["output"], par["output_compression"])
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], data, par["output_compression"]
)
logger.info("Finished")
