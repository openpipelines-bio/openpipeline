import scanpy as sc
import mudata as mu
import anndata as ad
import sys

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "base": None,
    "modality": "rna",
    "output_layer": "foo",
}
meta = {"name": "lognorm"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from input mudata", par["modality"])
data = mu.read_h5ad(par["input"], mod=par["modality"])
assert data.var_names.is_unique, "Expected var_names to be unique."

logger.info("Performing log transformation")
# Make our own copy with not a lot of data
# this avoid excessive memory usage and accidental overwrites
input_layer = data.layers[par["input_layer"]] if par["input_layer"] else data.X
data_for_scanpy = ad.AnnData(X=input_layer.copy())
sc.pp.log1p(
    data_for_scanpy,
    base=par["base"],
    layer=None,  # use X
    copy=False,
)  # allow overwrites in the copy that was made

# Scanpy will overwrite the input layer.
# So fetch input layer from the copy and use it to populate the output slot
if par["output_layer"]:
    data.layers[par["output_layer"]] = data_for_scanpy.X
else:
    data.X = data_for_scanpy.X
data.uns["log1p"] = data_for_scanpy.uns["log1p"].copy()

logger.info("Writing to file %s", par["output"])
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], data, par["output_compression"]
)
