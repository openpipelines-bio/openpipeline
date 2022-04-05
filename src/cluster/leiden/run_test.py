import subprocess
from os import path
import muon as mu

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tx_processing.h5mu"
output = "output.h5mu"

cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--obs_name", "fooleiden",
    "--output", output,
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists("output.h5mu"), "No output was created."

# read it with scanpy
data = mu.read_h5mu("output.h5mu")

# check whether gex was found
assert data.mod["rna"].var["feature_types"].unique() == [
    "Gene Expression"
], "Output X should only contain Gene Expression vars."

# check whether ab counts were found
assert "prot" in data.mod, 'Output should contain data.mod["rna"].'

# check whether leiden.custom.resolution was found
assert "fooleiden" in data.mod["rna"].obs.columns, 'Output should contain fooleiden.'
