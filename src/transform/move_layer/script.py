import sys
from mudata import read_h5ad
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

meta = {"resources_dir": "."}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s mudata from file %s", par["input"])
mod_data = read_h5ad(par["input"], mod=par["modality"])


logger.info(
    "Using input layer '%s'", "X" if not par["input_layer"] else par["input_layer"]
)
if par["input_layer"]:
    data_to_write = mod_data.layers[par["input_layer"]].copy()
    del mod_data.layers[par["input_layer"]]
else:
    data_to_write = mod_data.X
    mod_data.X = None

output_layer_setter = (
    partial(setattr, mod_data, "X")
    if not par["output_layer"]
    else partial(setitem, mod_data.layers, par["output_layer"])
)
output_layer_setter(data_to_write)

logger.info(
    "Writing output to file %s with compression %s",
    par["output"],
    par["output_compression"],
)

write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], mod_data, par["output_compression"]
)
