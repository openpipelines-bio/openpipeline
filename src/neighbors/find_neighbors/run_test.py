import subprocess
from os import path
import muon as mu
import logging
from sys import stdout

## VIASH START
meta = {
    'functionality_name': 'find_neighbors',
    'resources_dir': 'resources_test/'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
output = "output.h5mu"
cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--output", output,
    "--obsp_name_prefix", "foo"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if output was created.")
assert path.exists(output), "No output was created."

logger.info("> Reading mudata files.")
mu_input = mu.read_h5mu(input)
mu_output = mu.read_h5mu(output)

logger.info("> Check whether output contains right modalities.")
assert "rna" in mu_output.mod, 'Output should contain data.mod["prot"].'
assert "prot" in mu_output.mod, 'Output should contain data.mod["prot"].'

rna_in = mu_input.mod["rna"]
rna_out = mu_output.mod["rna"]
prot_in = mu_input.mod["prot"]
prot_out = mu_output.mod["prot"]

logger.info("> Check shape of outputs.")
assert rna_in.shape == rna_out.shape, "Should have same shape as before"
assert prot_in.shape == prot_out.shape, "Should have same shape as before"

logger.info("> Check existence of output fields.")
assert "foo_connectivities" in rna_out.obsp, "Output should have .obsp['foo_connectivities']"
assert "foo_distances" in rna_out.obsp, "Output should have .obsp['foo_distances']"
assert "foo_connectivities" not in rna_in.obsp, "Input should not have .obsp['foo_connectivities']"
assert "foo_distances" not in rna_in.obsp, "Input should not have .obsp['foo_distances']"