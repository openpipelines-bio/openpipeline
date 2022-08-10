import scanpy as sc
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
    "output": "foo.h5ad",
    "gex_only": False
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading %s.", par["input"])
data = sc.read_10x_h5(par["input"], gex_only=par["gex_only"])

logger.info("Making unique.")
data.var_names_make_unique()

logger.info("Writing to %s.", par["output"])
data.write_h5ad(par["output"])
logger.info("Finished")