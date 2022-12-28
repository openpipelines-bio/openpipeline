import os
import subprocess
import scanpy as sc
import mudata as mu
import logging
from sys import stdout
## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

orig_input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
input_path = "lognormed.h5mu"
output_path = "output.h5mu"

logger.info("> Prepare test data")
mu_in = mu.read_h5mu(orig_input_path)
rna_in = mu_in.mod["rna"]
assert "filter_with_hvg" not in rna_in.var.columns
sc.pp.log1p(rna_in)
mu_in.write_h5mu(input_path)


logger.info("> Run component")
out = subprocess.check_output(
    [
        f"./{meta['functionality_name']}",
        "--input", input_path,
        "--output", output_path
    ]
).decode("utf-8")

logger.info("> Check output file exists")
assert os.path.exists(output_path)
data = mu.read_h5mu(output_path)

logger.info("> Check whether column has been added")
assert "filter_with_hvg" in data.mod["rna"].var.columns