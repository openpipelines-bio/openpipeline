import subprocess
from os import path
import mudata as mu
import logging
from sys import stdout

## VIASH START
meta = {
    'resources_dir': 'resources_test'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("> Generate h5ad files")
mdata = mu.read_h5mu(meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu")
input_rna = "pbmc_1k_protein_v3_rna.h5ad"
mdata.mod["rna"].write_h5ad(input_rna)
input_prot = "pbmc_1k_protein_v3_prot.h5ad"
mdata.mod["prot"].write_h5ad(input_prot)

logger.info("> Run the command")
output = "output.h5mu"
cmd_pars = [
    meta["executable"],
    "--modality", "rna",
    "--input", input_rna,
    "--modality", "prot",
    "--input", input_prot,
    "--output", output,
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")
assert path.exists(output), "No output was created."

logger.info("> Read the output")
mdata2 = mu.read_h5mu(output)

assert "rna" in mdata2.mod, "Resulting mudata should contain rna modality"
assert "prot" in mdata2.mod, "Resulting mudata should contain rna modality"

logger.info("> All tests passed")