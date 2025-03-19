import mudata as mu
import sys

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "output_compression": "gzip",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# TODO: Merge modalities into one layer

logger.info("Reading input h5mu file %s, modality %s", par["input"], par["modality"])
adat = mu.read_h5ad(par["input"], mod=par["modality"])

logger.info("Writing to %s.", par["output"])
adat.write_h5ad(par["output"], compression=par["output_compression"])

logger.info("Finished")
