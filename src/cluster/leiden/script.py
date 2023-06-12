import logging
from sys import stdout
import mudata as mu
import scanpy as sc
from multiprocessing import Pool
from functools import partial

## VIASH START
par = {
    "input": "./resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "output_format": "h5mu",
    "obs_name_prefix": "leiden",
    "resolution": [1, 0.25],
    "obsp_connectivities": None,
    "uns_name_prefix": "leiden",
    "output_compression": "gzip"
}

meta = {
    'cpus': 4
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
    obs_key = f"{par['obs_name_prefix']}_{resolution}"
    adata_out = sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=obs_key,
        obsp=par['obsp_connectivities'],
        copy=True
        )
    return resolution, adata_out.obs[obs_key], adata_out.uns["leiden"]

logger.info("Processing modality '%s'.", par['modality'])
data = mdata.mod[par['modality']]
results = list(map(partial(run_single_resolution, data), par["resolution"]))
for resolution, result, result_parameters in results:
    data.obs[f"{par['obs_name_prefix']}_{resolution}"] = result
    data.uns[f"{par['uns_name_prefix']}_{resolution}"] = result_parameters
logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
logger.info("Finished.")
