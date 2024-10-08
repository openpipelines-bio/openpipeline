import sys
from mudata import read_h5mu

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "mod1",
    "input_obs_key": "Obs",
    "output_obs_key": "Obs_copy",
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

if not par["input_obs_key"] == par["output_obs_key"]:
    logger.info(f"Copying .obs key {par['input_obs_key']} to {par['output_obs_key']}")
    adata.obs[par["output_obs_key"]] = adata.obs[par["input_obs_key"]].copy()

logger.info("Write output to mudata file")
mdata.write_h5mu(par['output'], compression=par["output_compression"])
