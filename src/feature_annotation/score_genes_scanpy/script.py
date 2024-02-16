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

# START TEMPORARY WORKAROUND read_gene_list
# reason: resources aren't available when using Nextflow fusion

# from helper import read_gene_list
from typing import List, Dict, Any, Optional

def read_gene_list(
        par: Dict[str, Any],
        gene_names: List[str],
        list_key: str,
        file_key: str,
        required: bool = True) -> Optional[List[str]]:
    """
    Reads a gene list from the parameters and returns it as a list of strings.
    """

    # check whether one or the other was provided, if required
    if required and not par[list_key] and not par[file_key]:
        raise ValueError(f"Either --{list_key} or --{file_key} must be set")

    # read gene list from parameters
    list_of_genes = par[list_key] if par[list_key] else []

    # read gene list from file
    if par[file_key]:
        with open(par[file_key]) as file:
            file_genes = [x.strip() for x in file]
        list_of_genes.extend(file_genes)

    # check for missing genes
    if not par["allow_missing_genes"] and list_of_genes:
        missing = set(list_of_genes).difference(gene_names)
        if missing:
            raise ValueError(f"The follow genes are missing from the input dataset: {missing}")

    # return gene list
    if list_of_genes:
        return list_of_genes
    elif required:
        raise ValueError(f"No genes detected in --{list_key} or --{file_key}")
    else:
        return None

# END TEMPORARY WORKAROUND read_gene_list

# read data
mdata = mu.read(par["input"])
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
