import scanpy as sc
import mudata as mu
import logging
from sys import stdout
## VIASH START
par = {
    "input": "work/d9/3adbd080e0de618d44b59b1ec81685/run.output.h5mu",
    "output": "output.h5mu",
    "target_sum": 10000,
    "modality": "rna",
    "exclude_highly_expressed": False
}
meta = {"functionality_name": "lognorm"}
## VIASH END


logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

logger.info(par)

mod = par["modality"]
logger.info("Performing total normalization on modality %s", mod)
dat = mdata.mod[mod]
if par['input_layer'] and not par['input_layer'] in dat.layers.keys():
    raise ValueError(f"Input layer {par['input_layer']} not found in {mod}")
output_data = sc.pp.normalize_total(dat,
                                    layer=par["input_layer"],
                                    copy=True if par["output_layer"] else False)

if output_data:
    dat.layers[par["output_layer"]] = output_data.X

logger.info("Writing to file")
mdata.write_h5mu(filename=par["output"], compression="gzip")
