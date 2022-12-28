import mudata as mu
import anndata
import logging
from sys import stdout

## VIASH START
par = {
    "input": ["resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad"],
    "modality": ["rna"],
    "output": "output.h5mu",
    "compression": "gzip",
    "conversions_obsm": '{"counts_antibody":"prot", "counts_custom": "custom"}',
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

assert len(par["input"]) == len(par["modality"]), "Number of input files should be the same length as the number of modalities"

logger.info("Reading input files")
data = { key: anndata.read_h5ad(path) for key, path in zip(par["modality"], par["input"]) }

try:
    data.var_names_make_unique()
except:
    pass

logger.info("Converting to mudata")
mudata = mu.MuData(data)

try:
    mudata.var_names_make_unique()
except:
    pass

logger.info("Writing to %s.", par['output'])
mudata.write_h5mu(par["output"], compression="gzip")

logger.info("Finished")