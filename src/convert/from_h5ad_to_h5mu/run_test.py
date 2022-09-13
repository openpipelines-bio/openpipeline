import subprocess
from os import path
import muon as mu
import logging
from sys import stdout

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad"
output = "output.h5mu"

logger.info("> Run the command")
cmd_pars = [
    meta["executable"],
    "--input", input,
    "--output", output,
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")
assert path.exists(output), "No output was created."

logger.info("> Read the output")
data = mu.read_h5mu(output)

logger.info("> All tests passed")