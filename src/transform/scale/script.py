import sys
from mudata import read_h5ad
import scanpy
from functools import partial
from operator import setitem

## VIASH START
par = {
    "input": "../../../resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "max_value": None,
    "zero_center": True,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def main():
    logger.info(
        "Reading modality %s from .h5mu file: %s", par["modality"], par["input"]
    )
    data = read_h5ad(par["input"], mod=par["modality"])
    logger.info("Scaling modality %s.", par["modality"])
    scanpy_output = scanpy.pp.scale(
        data,
        layer=par["input_layer"],
        zero_center=par["zero_center"],
        max_value=par["max_value"],
        copy=True,
    )
    output_layer_setter = (
        partial(setattr, data, "X")
        if not par["output_layer"]
        else partial(setitem, data.layers, par["output_layer"])
    )
    output_layer_setter(
        scanpy_output.X
        if not par["input_layer"]
        else scanpy_output.layers[par["input_layer"]]
    )
    logger.info(
        "Writing to %s with compression %s", par["output"], par["output_compression"]
    )
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], data, par["output_compression"]
    )
    logger.info("Finished")


if __name__ == "__main__":
    main()
