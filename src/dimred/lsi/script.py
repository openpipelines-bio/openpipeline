import scanpy as sc
import muon as mu
import mudata as md
import logging
from anndata import AnnData
from sys import stdout
import numpy as np
import sys

## VIASH START
par = {
    "num_components": 50, # number of components to calculate with SVD
    "scale_embeddings": True, # scale embeddings to zero mean and unit variance
    "modality": "atac", # on which modality the LSI should be run 
    "layer": None, # on which layer to run the LSI, if None, will run it on anndata.X 
    "var_input": "highly_variable", # column in anndata.var of the highly variable features 
    #parameters for the output: 
    "overwrite": True, 
    "obsm_output": "X_lsi", # name of the slot in .obsm
    "varm_output": "LSI", #name of the slot in .varm
    "uns_output": "lsi", # name of the slot in .uns
    "output": "output.h5mu", #path to the output file 
    "output_compression": "gzip"
}
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


#1.read in mudata
logger.info("Reading %s.", par["input"])
mdata = md.read_h5mu(par["input"])

#2. subset on modality
if par["modality"] not in mdata.mod:
    raise ValueError(f" Modality '{par['modality']}' was not found in mudata {par['input']}.")
logger.info("Computing LSI for modality '%s'", par['modality'])
adata = mdata.mod[par['modality']]

#3. Specify layer
if par['layer'] and par["layer"] not in adata.layers:
    raise ValueError(f" Layer '{par['layer']}' was not found in modality '{par['modality']}'.")
layer = adata.X if not par['layer'] else adata.layers[par['layer']]
adata_input_layer = AnnData(layer)
adata_input_layer.var.index = adata.var.index

#4. Subset on highly variable features if applicable
use_highly_variable = False
if par["var_input"]:
    if not par["var_input"] in adata.var.columns:
        raise ValueError(f"Requested to use .var column '{par['var_input']}' "
                         "as a selection of genes to run the LSI on, "
                         f"but the column is not available for modality '{par['modality']}'")
    use_highly_variable = True
    #subset anndata on highly variable 
    adata_input_layer = adata_input_layer[:, adata.var[par["var_input"]]] 

#5. Run LSI
mu.atac.tl.lsi(adata_input_layer, scale_embeddings = par["scale_embeddings"], n_comps = par["num_components"])



#6. Store output in object
check_exist_dict = {
    "obsm_output": ("obs"),
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
if use_highly_variable: 
    adata.varm[par["varm_output"]] = np.zeros(shape=(adata.n_vars, adata_input_layer.varm["LSI"].shape[1]))
    adata.varm[par["varm_output"]][adata.var[par["var_input"]]] = adata_input_layer.varm['LSI']
else:
    adata.varm[par["varm_output"]] = adata_input_layer.varm['LSI']

logger.info("Writing to %s.", par["output"])
mdata.write(filename = par["output"], compression=par["output_compression"])

logger.info("Finished")
