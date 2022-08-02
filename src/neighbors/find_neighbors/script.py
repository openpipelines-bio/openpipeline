import muon as mu
import scanpy as sc
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "metric": 'cosine',
    "num_neighbors": 15,
    "modality": ["rna"],
    "obsp_name_prefix": "foo"
}
meta = {"functionality_name": "lognorm"}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

for mod in par["modality"]:
    logger.info("Computing a neighborhood graph on modality %s", mod)
    sc.pp.neighbors(
        mdata.mod[mod],
        n_neighbors=par["num_neighbors"], 
        metric=par["metric"],
        key_added=par["obsp_name_prefix"]
    )

logger.info("Writing to %s", par["output"])
mdata.write(filename=par["output"])
