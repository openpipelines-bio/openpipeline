import subprocess
from os import path
import scanpy as sc

cmd_pars = ["./convert_10x_to_h5ad", "--input", "pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",  "--output=output.h5ad"]
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists("output.h5ad"), "No output was created."

# read it with scanpy
data = sc.read_h5ad("output.h5ad")

# check whether gex was found
assert data.var["feature_types"].unique() == ["Gene Expression"], "Output X should only contain Gene Expression vars."

# check whether ab counts were found
assert "counts_antibody" in data.obsm_keys(), "Output should contain adata.obsm[\"counts_antibody\"]."

# check whether gene was found
assert "CD3_TotalSeqB" in data.obsm["counts_antibody"].columns, "Output should contain antibody column \"CD3_TotalSeqB\"."

