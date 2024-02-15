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
    "gene_list": ["UBE2C", "BIRC5", "TPX2"],
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
meta = {
    "resources_dir": "src/feature_annotation/score_genes_scanpy"
}
## VIASH END

sys.path.append(meta["resources_dir"])

from helper import read_gene_list

# read data
mdata = mu.read(f"{par["input"]}")
input_adata = mdata.mod[par["modality"]]

if par["var_gene_names"]:
    gene_names = list(input_adata.var[par["var_gene_names"]])
else:
    gene_names = list(input_adata.var.index)

# check if var index is unique
if not input_adata.var.index.is_unique:
    raise ValueError("var index is not unique")


# read gene list
gene_list = read_gene_list(
    par,
    gene_names,
    "gene_list",
    "gene_list_file"
)
gene_pool = read_gene_list(
    par,
    gene_names,
    "gene_pool",
    "gene_pool_file",
    required=False
)

# find matching index names for given genes
gene_list_index = input_adata.var.index[[gene in gene_list for gene in gene_names]]
gene_pool_index = input_adata.var.index[[gene in gene_pool for gene in gene_names]] if gene_pool else None

# create input data for scanpy
layer_data = input_adata.layers[par["input_layer"]] if par["input_layer"] else input_adata.X
adata_scanpy = ad.AnnData(
    X=layer_data,
    obs=pd.DataFrame(index=input_adata.obs.index),
    var=pd.DataFrame(index=input_adata.var.index)
)

# run score_genes
sc.tl.score_genes(
    adata_scanpy,
    gene_list=gene_list_index,
    gene_pool=gene_pool_index,
    ctrl_size=par["ctrl_size"],
    n_bins=par["n_bins"],
    random_state=par["random_state"]
)

# copy results to mudata
input_adata.obs[par["obs_score"]] = adata_scanpy.obs["score"]

# write output to mudata
mdata.write(par["output"], compression=par["output_compression"])
