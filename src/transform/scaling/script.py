import sys
from mudata import read_h5mu
import scanpy

## VIASH START
par = {
    "input": "../../../resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "max_value": None,
    "zero_center": True
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

def main():
    logger.info(f'Reading .h5mu file: {par["input"]}')
    mudata = read_h5mu(par["input"])
    mod = par["modality"]
    data = mudata.mod[mod]

    logger.info("Scaling modality: %s", mod)
    scanpy.pp.scale(data,
                    zero_center=par["zero_center"],
                    max_value=par["max_value"])

    logger.info("Writing to %s", par["output"])
    mudata.write_h5mu(filename=par["output"], compression=par["output_compression"])
    logger.info("Finished")

if __name__ == "__main__":
    main()