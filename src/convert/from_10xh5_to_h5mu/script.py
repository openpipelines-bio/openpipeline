import muon as mu
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.h5",
    "output": "5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5mu",
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading %s.", par["input"])
mdata = mu.read_10x_h5(par["input"])

logger.info("Renaming keys.")
for adata in mdata.mod.items:
    adata.var.rename(columns={'gene_ids': 'gene_id', 'feature_types': 'feature_type'}, inplace=True)

logger.info("Making unique.")
mdata.var_names_make_unique()

logger.info("Writing %s.", par["output"])
mdata.write_h5mu(par["output"])
