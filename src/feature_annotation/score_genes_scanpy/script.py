import scanpy as sc
import anndata as ad
import pandas as pd
import mudata as mu
import sys

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_layer": "log_normalized",
    "gene_list_file": None,
    "gene_list": ["MCM5", "PCNA", "TYMS"],
    "gene_pool_file": None,
    "gene_pool": None,
    "var_gene_names": "gene_symbol",
    "output": "output.h5mu",
    "ctrl_size": 50,
    "n_bins": 25,
    "obs_score": "score",
    "random_state": 0,
    "output_compression": "gzip",
    "allow_missing_genes": False,
}
meta = {"resources_dir": "src/feature_annotation/score_genes_scanpy"}
## VIASH END

sys.path.append(meta["resources_dir"])
from helper import read_gene_list

# read data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]

gene_names_index = (
    input_adata.var[par["var_gene_names"]]
    if par["var_gene_names"]
    else input_adata.var_names
)
gene_names = pd.Series(input_adata.var_names, index=gene_names_index)

# check if var index is unique
# input.var[par["var_gene_names"]] is mapped to var index, but may not contain unique values
if not input_adata.var.index.is_unique:
    raise ValueError("var index is not unique")

# read gene list
gene_list = read_gene_list(par, gene_names.index, "gene_list", "gene_list_file")
gene_pool = read_gene_list(
    par, gene_names.index, "gene_pool", "gene_pool_file", required=False
)

# find matching index names for given genes
gene_list_index = gene_names.loc[gene_list].tolist()
gene_pool_index = gene_names.loc[gene_pool].tolist() if gene_pool else None

# create input data for scanpy
if par["input_layer"]:
    layer_data = input_adata.layers[par["input_layer"]].copy()
else:
    layer_data = input_adata.X.copy()
adata_scanpy = ad.AnnData(
    X=layer_data,
    obs=pd.DataFrame(index=input_adata.obs.index),
    var=pd.DataFrame(index=input_adata.var.index),
)

# run score_genes
sc.tl.score_genes(
    adata_scanpy,
    gene_list=gene_list_index,
    gene_pool=gene_pool_index,
    ctrl_size=par["ctrl_size"],
    n_bins=par["n_bins"],
    random_state=par["random_state"],
)

# copy results to mudata
assert all(adata_scanpy.obs.index == input_adata.obs.index), (
    "index mismatch between input adata and scanpy output adata"
)
input_adata.obs[par["obs_score"]] = adata_scanpy.obs["score"]

# write output to mudata
mdata.write(par["output"], compression=par["output_compression"])
