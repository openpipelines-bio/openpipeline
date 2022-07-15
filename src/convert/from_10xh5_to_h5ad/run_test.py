import subprocess
from os import path
import scanpy as sc
import numpy as np

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5"

############### setting 1
output = "output1.h5ad"

cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--output", output,
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists(output), "No output was created."

# read it with scanpy
data = sc.read_h5ad(output)

# check whether gex was found
assert np.array_equal(data.var["feature_types"].unique(), ["Gene Expression", "Antibody Capture"]), "Output X should contain Gene Expression and ADT vars."

# check whether gene was found
assert "CD3_TotalSeqB" in data.var_names, 'Output should contain antibody column "CD3_TotalSeqB".'

############### setting 2
output = "output2.h5ad"

cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--output", output,
    "--gex_only"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists(output), "No output was created."

# read it with scanpy
data = sc.read_h5ad(output)

# check whether gex was found
assert np.array_equal(data.var["feature_types"].unique(), ["Gene Expression"]), "Output X should only contain Gene Expression vars."

# check whether gene was found
assert "OR4F5" in data.var_names, 'Output should contain antibody column "CD3_TotalSeqB".'
