import scanpy as sc
import mudata as mu
import anndata as ad
import pandas as pd
import sys

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_layer": "log_normalized",
    "var_gene_names": "gene_symbol",
    "s_genes": ["MCM5", "PCNA", "TYMS"],
    "s_genes_file": None,
    "g2m_genes": ["UBE2C", "BIRC5", "TPX2"],
    "g2m_genes_file": None,
    "gene_pool": [],
    "gene_pool_file": None,
    "output": "output.h5mu",
    "obs_phase": "phase",
    "obs_s_score": "S_score",
    "obs_g2m_score": "G2M_score",
    "n_bins": 25,
    "random_state": 0,
    "output_compression": "gzip",
    "allow_missing_genes": False
}
meta = {
    "resources_dir": "src/feature_annotation/score_genes_scanpy"
}
## VIASH END

# import helper functions
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
# input.var[par["var_gene_names"]] is mapped to var index, but may not contain unique values
if not input_adata.var.index.is_unique:
    raise ValueError("var index is not unique")

# read gene lists
s_genes = read_gene_list(par, gene_names, "s_genes", "s_genes_file")
g2m_genes = read_gene_list(par, gene_names, "g2m_genes", "g2m_genes_file")
gene_pool = read_gene_list(par, gene_names, "gene_pool", "gene_pool_file", required=False)

# find matching index names for given genes
gene_names_index = input_adata.var[par["var_gene_names"]] if par["var_gene_names"] else input_adata.var_names
gene_names = pd.Series(input_adata.var_names, index=gene_names_index)

g2m_index = gene_names.loc[g2m_genes].tolist()
s_index = gene_names.loc[s_genes].tolist()
gene_pool_index = gene_names.loc[gene_pool].tolist() if gene_pool else None

# create input data for scanpy
if par["input_layer"]:
    X_data = input_adata.layers[par["input_layer"]].copy()
else:
    X_data = input_adata.X.copy()
adata_scanpy = ad.AnnData(
    X=X_data,
    obs=pd.DataFrame(index=input_adata.obs.index),
    var=pd.DataFrame(index=input_adata.var.index)
)

# run score_genes_cell_cycle
sc.tl.score_genes_cell_cycle(
    adata_scanpy,
    s_genes=s_index,
    g2m_genes=g2m_index,
    gene_pool=gene_pool_index,
    n_bins=par["n_bins"],
    random_state=par["random_state"]
)

# copy results to mudata
output_slot_mapping = {
    par["obs_s_score"]: "S_score",
    par["obs_g2m_score"]: "G2M_score",
    par["obs_phase"]: "phase"
}
assert all(adata_scanpy.obs.index == input_adata.obs.index), "index mismatch between input adata and scanpy output adata"
for dest, orig in output_slot_mapping.items():
    input_adata.obs[dest] = adata_scanpy.obs[orig]

# write output to mudata
mdata.write(par["output"], compression=par["output_compression"])
