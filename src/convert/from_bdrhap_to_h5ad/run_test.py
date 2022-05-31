import subprocess
from os import path
import anndata as ad
import numpy as np

input = meta["resources_dir"] + "/bd_rhapsody_wta_test/processed/bdrhap_out"
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
data = ad.read_h5ad(output)

# check whether gex was found
assert np.array_equal(data.var["feature_types"].unique(), ["Gene Expression"]), "Output X should only contain Gene Expression vars."

# check whether gene was found
assert "PDE4DIP" in data.var_names, 'Output should contain gex column "PDE4DIP".'
