import muon as mu
import mudata as md
from anndata import AnnData
import numpy as np
import sys


## VIASH START
par = {
    "num_components": 50, # number of components to calculate with SVD
    "scale_embeddings": True, # scale embeddings to zero mean and unit variance
    "modality": "atac", # on which modality the LSI should be run 
    "layer": None, # on which layer to run the LSI, if None, will run it on anndata.X 
    "var_input": None, # column in anndata.var of the highly variable features

    "overwrite": True, 
    "obsm_output": "X_lsi",
    "varm_output": "LSI",
    "uns_output": "lsi",
    "output": "output.h5mu",
    "output_compression": "gzip"
}
## VIASH END


sys.path.append(meta["resources_dir"])
from subset_vars import subset_vars


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


#1.read in mudata
logger.info("Reading %s.", par["input"])
mdata = md.read_h5mu(par["input"])

#2. subset on modality
if par["modality"] not in mdata.mod:
    raise ValueError(f"Modality '{par['modality']}' was not found in mudata {par['input']}.")
adata = mdata.mod[par['modality']]


#3. Specify layer
if par['layer'] and par["layer"] not in adata.layers:
    raise ValueError(f"Layer '{par['layer']}' was not found in modality '{par['modality']}'.")
layer = adata.X if not par['layer'] else adata.layers[par['layer']]
adata_input_layer = AnnData(layer, var=adata.var)


if not par["layer"]:
    logger.info("Using modality '%s' and adata.X for LSI computation", par['modality'])
else:
    logger.info("Using modality '%s' and layer '%s' for LSI computation", par['modality'], par["layer"])


#4. Subset on highly variable features if applicable
if par["var_input"]:
    adata_input_layer = subset_vars(adata_input_layer, par["var_input"])



#5. Run LSI
logger.info("Computing %s LSI components on %s features", par["num_components"], adata_input_layer.X.shape[1])
mu.atac.tl.lsi(adata_input_layer, scale_embeddings = par["scale_embeddings"], n_comps = par["num_components"])



#6. Store output in object
check_exist_dict = {
    "obsm_output": ("obsm"),
    "varm_output": ("varm"),
    "uns_output": ("uns")
}
for parameter_name, field in check_exist_dict.items():
    if par[parameter_name] in getattr(adata, field):
        if not par["overwrite"]:
            raise ValueError(f"Requested to create field {par[parameter_name]} in .{field} "
                            f"for modality {par['modality']}, but field already exists.")
        del getattr(adata, field)[par[parameter_name]]

adata.obsm[par["obsm_output"]] = adata_input_layer.obsm['X_lsi']
adata.uns[par["uns_output"]] = adata_input_layer.uns['lsi']
if par["var_input"]:
    adata.varm[par["varm_output"]] = np.zeros(shape=(adata.n_vars, adata_input_layer.varm["LSI"].shape[1]))
    adata.varm[par["varm_output"]][adata.var[par["var_input"]]] = adata_input_layer.varm['LSI']
else:
    adata.varm[par["varm_output"]] = adata_input_layer.varm['LSI']

logger.info("Writing to %s.", par["output"])
mdata.write(filename = par["output"], compression=par["output_compression"])

logger.info("Finished")
