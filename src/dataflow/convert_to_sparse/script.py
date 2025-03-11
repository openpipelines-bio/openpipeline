import mudata as mu
import sys
from scipy.sparse import issparse, csr_matrix

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "layer": None,
    "output": "output.h5mu",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/utils"}

mdata = mu.read_h5mu(par["input"])
mdata.mod["rna"].X = mdata.mod["rna"].X.toarray()
mdata.write_h5mu("dense_input.h5mu")
par["input"] = "dense_input.h5mu"
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# TODO: Merge modalities into one layer

logger.info("Reading input h5mu file %s, modality %s", par["input"], par["modality"])
mdata = mu.read_h5mu(par["input"])
adat = mdata.mod[par["modality"]].copy()

layer = adat.layer[par["layer"]] if par["layer"] else adat.X
layer_name = par["layer"] if par["layer"] else ".X"

if issparse(layer):
    logger.warning(f"{layer_name} is already sparse, skipping conversion.")

else:
    logger.info(f"Converting {layer_name} to sparse matrix.")
    layer_sparse = csr_matrix(layer)
    if par["layer"]:
        adat.layer[par["layer"]] = layer_sparse
    else:
        adat.X = layer_sparse

    mdata.mod[par["modality"]] = adat

# Write mudata output
logger.info("Writing output data")
mdata.write_h5mu(par["output"], compression=par["output_compression"])

logger.info("Finished")
