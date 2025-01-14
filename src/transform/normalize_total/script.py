import sys
import scanpy as sc
import mudata as mu

## VIASH START
par = {
    "input": "work/d9/3adbd080e0de618d44b59b1ec81685/run.output.h5mu",
    "output": "output.h5mu",
    "target_sum": 10000,
    "modality": "rna",
    "exclude_highly_expressed": False,
}
meta = {"name": "lognorm"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from %s", par["modality"], par["input"])
dat = mu.read_h5ad(par["input"], mod=par["modality"])
assert dat.var_names.is_unique, "The var_names of the input modality must be be unique."

logger.info(par)

logger.info("Performing total normalization.")
if par["input_layer"] and par["input_layer"] not in dat.layers.keys():
    raise ValueError(f"Input layer {par['input_layer']} not found in {par['modality']}")
output_data = sc.pp.normalize_total(
    dat,
    layer=par["input_layer"],
    target_sum=par["target_sum"],
    copy=True if par["output_layer"] else False,
)

if output_data:
    result = (
        output_data.X
        if not par["input_layer"]
        else output_data.layers[par["input_layer"]]
    )
    dat.layers[par["output_layer"]] = result

logger.info(
    "Writing to file to %s with compression %s",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], dat, par["output_compression"]
)
