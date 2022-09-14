import logging
from sys import stdout
import muon as mu
import scanpy as sc

## VIASH START
par = {
    "input": "work/c3/a5a0813d70d0141193748e3baf9a58/pbmc_1k_protein_v3_mms.find_neighbors.output.h5mu",
    "output": "output.h5mu",
    "modality": ["rna"],
    "output_format": "h5mu",
    "obs_name": "leiden",
    "resolution": 0.25,
    "obsp_connectivities": "connectivities"
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading %s.", par["input"])
mdata = mu.read_h5mu(par["input"])

for mod in par['modality']:
    logger.info("Processing modality '%s'.", mod)
    data = mdata.mod[mod]
    sc.tl.leiden(
        data,
        resolution=par["resolution"],
        key_added=par["obs_name"],
        obsp=par["obsp_connectivities"]
    )

logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"])
logger.info("Finished.")
