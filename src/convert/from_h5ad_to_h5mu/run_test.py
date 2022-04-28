import subprocess
from os import path
import muon as mu

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad"
output = "output.h5mu"

print("> Run the command")
cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--output", output,
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

print("> Check if file exists")
assert path.exists(output), "No output was created."

print("> Read the output")
data = mu.read_h5mu(output)

print("> All tests passed")