import sys
from mudata import read_h5mu

## VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "modality": "rna",
    "input_var_key": None,
    "output_var_key": "index_copy",
    "output": "output.h5mu",
    "output_compression": "gzip",
    "disable_raise_on_identical_keys": False
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


def duplicate_var(adata, input_key, output_key):
    if input_key:
        logger.info(f"Copying .var key {input_key} to {output_key}")
        adata.var[output_key] = adata.var[input_key].copy()
    else:
        logger.info(f"Copying .var index to {output_key}")
        adata.var[output_key] = adata.var.index.copy()


if not par["output_var_key"] in adata.var:
    duplicate_var(adata, par["input_var_key"], par["output_var_key"])

else:
    if par["overwrite_existing_key"]:
        logger.warning(f"--output_var_key already exists: `{par['output_var_key']}`. Data in `{par['output_var_key']}` .var column will be overwritten.")
        duplicate_var(adata, par["input_var_key"], par["output_var_key"])
    else:
        raise ValueError(f"--output_var_key already exists: `{par['output_var_key']}`. Data can not be duplicated.")

logger.info("Write output to mudata file")

mdata.write_h5mu(par['output'], compression=par["output_compression"])
