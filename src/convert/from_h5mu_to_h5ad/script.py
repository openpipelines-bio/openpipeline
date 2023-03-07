import mudata as mu
import logging
from sys import stdout

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "output_compression": "gzip",
}
## VIASH END

# TODO: Merge modalities into one layer

logger.info("Reading input h5mu file")
dat = mu.read_h5mu(par["input"])

logger.info("Converting to h5ad")
adat = dat.mod[par["modality"]]

logger.info("Writing to %s.", par['output'])
adat.write_h5ad(par["output"], compression=par["output_compression"])

logger.info("Finished")