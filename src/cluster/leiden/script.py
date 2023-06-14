import logging
from sys import stdout
import mudata as mu
import pandas as pd
import scanpy as sc
from functools import partial

## VIASH START
par = {
    "input": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "output_format": "h5mu",
    "obsm_name": "leiden",
    "resolution": [1, 0.25],
    "obsp_connectivities": None,
    "uns_name": "leiden",
    "output_compression": "gzip"
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


def run_single_resolution(adata, resolution):
    adata_out = sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=str(resolution),
        obsp=par['obsp_connectivities'],
        copy=True
        )
    return adata_out.obs[str(resolution)]

logger.info("Processing modality '%s'.", par['modality'])
data = mdata.mod[par['modality']]
results = {str(resolution): run_single_resolution(data, resolution) for resolution in par["resolution"]}
data.obsm[par["obsm_name"]] = pd.DataFrame(results)
logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
logger.info("Finished.")
