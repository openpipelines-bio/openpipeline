import mudata as mu
import scanpy as sc
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix",
    "output": "foo.h5mu",
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading %s.", par["input"])
adata = sc.read_10x_mtx(par["input"], gex_only=False)

logger.info("Renaming keys.")
adata.var = adata.var\
  .rename_axis("gene_symbol")\
  .reset_index()\
  .set_index("gene_ids")

# generate output
logger.info("Convert to mudata")
mdata = mu.MuData(adata)

# override root .obs
mdata.obs = adata.obs

# write output
logger.info("Writing %s", par["output"])
mdata.write_h5mu(par["output"], compression=par["output_compression"])