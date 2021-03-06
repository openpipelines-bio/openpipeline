import subprocess
from os import path
import muon as mu

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

input = meta["resources_dir"] + "pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
output = "output.h5mu"

cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--output", output,
    "--output_key", "foo",
    "--num_components", "26"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists("output.h5mu"), "No output was created."

# read it with scanpy
data = mu.read_h5mu("output.h5mu")

# check whether pca was found
assert "X_foo" in data.mod["rna"].obsm, "Check whether output was found in .obsm"
assert data.mod["rna"].obsm["X_foo"].shape == (data.n_obs, 26), "Check shapes"