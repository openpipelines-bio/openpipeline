import mudata as mu
import anndata
import sys

## VIASH START
par = {
    "input": [
        "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad"
    ],
    "modality": ["rna"],
    "output": "output.h5mu",
    "output_compression": "gzip",
    "conversions_obsm": '{"counts_antibody":"prot", "counts_custom": "custom"}',
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

assert len(par["input"]) == len(
    par["modality"]
), "Number of input files should be the same length as the number of modalities"

logger.info("Reading input files")
data = {
    key: anndata.read_h5ad(path) for key, path in zip(par["modality"], par["input"])
}

logger.info("Converting to mudata")
mudata = mu.MuData(data)

try:
    mudata.var_names_make_unique()
except (TypeError, ValueError):
    pass

logger.info("Writing to %s.", par["output"])
mudata.write_h5mu(par["output"], compression=par["output_compression"])

logger.info("Finished")
