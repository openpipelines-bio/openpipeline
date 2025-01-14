import sys
from scanpy.external.pp import scanorama_integrate
from mudata import read_h5ad


### VIASH START
par = {}
meta = {"resources_dir": "src/utils/"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from input %s", par["modality"], par["input"])
mod = read_h5ad(par["input"], mod=par["modality"])

arguments = {
    "key": par["obs_batch"],
    "basis": par["obsm_input"],
    "adjusted_basis": par["obsm_output"],
    "knn": par["knn"],
    "alpha": par["alpha"],
    "sigma": par["sigma"],
    "approx": par["approx"],
    "batch_size": par["batch_size"],
}

logger.info(
    "Running scanorama with parameters: \n",
    "\n\t".join(
        [
            f"{argument_name}: {argument_value}"
            for argument_name, argument_value in arguments.items()
        ]
    ),
)
# Integration.
scanorama_integrate(mod, **arguments)
logger.info("Writing output to %s", par["output"])
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], mod, par["output_compression"]
)

logger.info("Finished!")
