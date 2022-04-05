import subprocess
from os import path
import muon as mu


## VIASH START
meta = {
    'functionality_name': 'find_neighbors',
    'resources_dir': 'resources_test/'
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tx_processing.h5mu"
output = "output.h5mu"
cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--output", output,
    "--obsp_name_prefix", "foo"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

print("> Check if output was created.")
assert path.exists(output), "No output was created."

print("> Reading mudata files.")
mu_input = mu.read_h5mu(input)
mu_output = mu.read_h5mu(output)

print("> Check whether output contains right modalities.")
assert "rna" in mu_output.mod, 'Output should contain data.mod["prot"].'
assert "prot" in mu_output.mod, 'Output should contain data.mod["prot"].'

rna_in = mu_input.mod["rna"]
rna_out = mu_output.mod["rna"]
prot_in = mu_input.mod["prot"]
prot_out = mu_output.mod["prot"]

print("> Check shape of outputs.")
assert rna_in.shape == rna_out.shape, "Should have same shape as before"
assert prot_in.shape == prot_out.shape, "Should have same shape as before"

print("> Check existence of output fields.")
assert "foo_connectivities" in rna_out.obsp, "Output should have .obsp['foo_connectivities']"
assert "foo_distances" in rna_out.obsp, "Output should have .obsp['foo_distances']"
assert "foo_connectivities" not in rna_in.obsp, "Input should not have .obsp['foo_connectivities']"
assert "foo_distances" not in rna_in.obsp, "Input should not have .obsp['foo_distances']"