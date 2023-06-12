import subprocess
from os import path
import mudata as mu

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
output = "output.h5mu"

cmd_pars = [
    meta["executable"],
    "--input", input,
    "--obs_name_prefix", "fooleiden",
    "--resolution", "1",
    "--resolution", "0.25",
    "--output", output,
    "--output_compression", "gzip",
    "--uns_name_prefix", "barleiden",
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists("output.h5mu"), "No output was created."

# read it with scanpy
data = mu.read_h5mu("output.h5mu")

# check whether leiden.custom.resolution was found
assert "fooleiden_1.0" in data.mod["rna"].obs.columns, 'Output should contain fooleiden_1.'
assert "fooleiden_0.25" in data.mod["rna"].obs.columns, 'Output should contain fooleiden_0.25.'

# check whether .uns can be found was found
assert "barleiden_1.0" in data.mod["rna"].uns, '.uns should contain barleiden_1.'
assert "barleiden_0.25" in data.mod["rna"].uns, '.uns should contain barleiden_0.25.'
