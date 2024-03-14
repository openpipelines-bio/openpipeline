import sys
from mudata import read_h5mu
from functools import partial
from operator import setitem

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "foo.h5mu",
    "modality": "rna",
    "input_layer": None,
    "output_layer": "test",
    "output_compression": None,
}

meta = {
    "resources_dir": "."
}
### VIASH END

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
input_file, modality = par["input"], par["modality"]
mdata = read_h5mu(input_file)
mod_data = mdata.mod[modality]


logger.info("Using input layer '%s'", "X" if not par["input_layer"] else par["input_layer"])
if par["input_layer"]:
    data_to_write = mod_data.layers[par["input_layer"]].copy()
    del mod_data.layers[par["input_layer"]]
else:
    data_to_write = mod_data.X
    mod_data.X = None

output_layer_setter = partial(setattr, mod_data, "X") \
                        if not par["output_layer"] \
                        else partial(setitem, mod_data.layers, par["output_layer"])
output_layer_setter(data_to_write)

logger.info("Write output to mudata file")
mdata.write_h5mu(par['output'], compression=par["output_compression"])

        

