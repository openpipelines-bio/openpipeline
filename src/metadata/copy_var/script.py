import sys
from mudata import read_h5mu

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "mod1",
    "input_var_key": None,
    "output_var_key": "Feat_copy",
    "output": "output.h5mu",
    "output_compression": "gzip"
}
meta = {
    "resources_dir": "src/metadata/copy_var"
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

if not par["input_var_key"] == par["output_var_key"]:
    if par["input_var_key"]:
        logger.info(f"Copying .var key {par['input_var_key']} to {par['output_var_key']}")
        adata.var[par["output_var_key"]] = adata.var[par["input_var_key"]].copy()
    else:
        logger.info(f"Copying .var index to {par['output_var_key']}")
        adata.var[par["output_var_key"]] = adata.var.index.copy()
logger.info("Write output to mudata file")

mdata.write_h5mu(par['output'], compression=par["output_compression"])
