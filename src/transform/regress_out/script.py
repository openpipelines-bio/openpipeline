import scanpy as sc
import mudata as mu
import multiprocessing
import sys

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obs_keys": [],
}
meta = {"functionality_name": "lognorm"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

if (
    par["obs_keys"] is not None
    and len(par["obs_keys"]) > 0
):
    mod = par["modality"]
    logger.info("Regress out variables on modality %s", mod)
    data = mdata.mod[mod]

    sc.pp.regress_out(
        data,
        keys=par["obs_keys"],
        n_jobs=multiprocessing.cpu_count() - 1
    )

logger.info("Writing to file")
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
