import subprocess
from os import path
import muon as mu

## VIASH START
meta = {
  'executable': 'target/docker/correction/cellbender_remove_background/cellbender_remove_background',
  'resources_dir': 'resources_test/pbmc_1k_protein_v3'
}
## VIASH END

file_raw = meta["resources_dir"] + "/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5"
file_input = "input.h5mu"
file_output = "output.h5mu"

print("> Check whether cellbender works when it should be working")
# read 10x h5 file and write as h5mu
mdat = mu.read_10x_h5(file_raw)
mdat = mdat[0:100000,] # subsample to reduce computational time
mdat.write_h5mu(file_input)

# run cellbender
cmd_pars = [
    meta["executable"],
    "--input", file_input,
    "--output", file_output,
    "--epochs", "5",
    "--output_compression", "gzip"
]
# todo: if cuda is available, add --cuda
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