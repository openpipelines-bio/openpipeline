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
meta = {"functionality_name": "lognorm"}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

mod = par["modality"]
logger.info("Performing log transformation on modality %s", mod)
data = mdata.mod[mod]

# Make our own copy with not a lot of data
# this avoid excessive memory usage and accidental overwrites 
input_layer = data.layers[par["input_layer"]] \
              if par["input_layer"] else data.X
data_for_scanpy = ad.AnnData(X=input_layer.copy())
sc.pp.log1p(data_for_scanpy,
            base=par["base"],
            layer=None, # use X
            copy=False) # allow overwrites in the copy that was made

# Scanpy will overwrite the input layer.
# So fetch input layer from the copy and use it to populate the output slot
if par["output_layer"]:
    data.layers[par["output_layer"]] = data_for_scanpy.X
else:
    data.X = data_for_scanpy.X
data.uns['log1p'] = data_for_scanpy.uns['log1p'].copy()

logger.info("Writing to file %s", par["output"])
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
