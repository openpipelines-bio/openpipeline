import sys
from mudata import read_h5mu

## VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "modality": "rna",
    "overwrite_existing_key": False,
    "input_obs_key": None,
    "output_obs_key": "index_copy",
    "output": "output.h5mu",
    "output_compression": "gzip"
}
meta = {
    "resources_dir": "src/metadata/copy_obs"
}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

logger.info("Read mudata from file")
mdata = read_h5mu(par['input'])
adata = mdata.mod[par['modality']]


def duplicate_obs(adata, input_key, output_key):
    if input_key:
        logger.info(f"Copying .obs key {input_key} to {output_key}")
        adata.obs[output_key] = adata.obs[input_key].copy()
    else:
        logger.info(f"Copying .obs index to {output_key}")
        adata.obs[output_key] = adata.obs.index.copy()

     
if not par["output_obs_key"] in adata.obs:
    duplicate_obs(adata, par["input_obs_key"], par["output_obs_key"])

else:
    if not par["overwrite_existing_key"]:
        raise ValueError(f"--output_obs_key already exists: `{par['output_obs_key']}`. Data can not be duplicated.")

    logger.warning(f"--output_obs_key already exists: `{par['output_obs_key']}`. Data in par['output_obs_key'] will be overwritten.")
    duplicate_obs(adata, par["input_obs_key"], par["output_obs_key"])


logger.info("Write output to mudata file")
mdata.write_h5mu(par['output'], compression=par["output_compression"])
