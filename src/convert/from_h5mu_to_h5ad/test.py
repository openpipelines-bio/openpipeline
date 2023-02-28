import subprocess
from os import path
from sys import stdout
import logging
import anndata as ad
import mudata as mu

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
meta = {
    'executable': 'target/docker/convert/from_h5mu_to_h5ad/from_h5mu_to_h5ad',
    'resources_dir': 'resources_test'
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"
output = "output.h5ad"

logger.info("> Run the command")
output = "output.h5mu"
cmd_pars = [
    meta["executable"],
    "--modality", "rna",
    "--input", input,
    "--output", output,
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")
assert path.exists(output), "No output was created."

logger.info("> Read the output")
adata = ad.read_h5ad(output)

logger.info("> Compare against expected output")
mdata = mu.read_h5mu(input)
assert adata.n_obs == mdata.mod["rna"].n_obs
assert adata.n_vars == mdata.mod["rna"].n_vars

logger.info("> All tests passed")