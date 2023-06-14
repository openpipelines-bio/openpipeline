import subprocess
from os import path
import mudata as mu
import numpy as np
import logging
from sys import stdout

## VIASH START
meta = {
    'functionality_name': 'lognorm',
    'resources_dir': 'resources_test/'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

orig_input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
input = "input.h5mu"
output = "output.h5mu"

# add column to obs
mu_orig_input = mu.read_h5mu(orig_input)
mu_orig_input.obs["var"] = np.random.rand(mu_orig_input.n_obs)
mu_orig_input.mod["rna"].obs["var"] = mu_orig_input.obs["var"]
mu_orig_input.mod["prot"].obs["var"] = mu_orig_input.obs["var"]
mu_orig_input.write_h5mu(input)

# execute command
cmd_pars = [
    f"./{meta['functionality_name']}",
    "--input", input,
    "--output", output,
    "--obs_keys", "var",
    "--output_compression", "gzip"
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

logger.info("> Check if expression has changed.")
assert np.mean(rna_in.X) != np.mean(rna_out.X), "Expression should have changed"

logger.info(">> All tests succeeded!")