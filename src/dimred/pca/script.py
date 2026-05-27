import scanpy as sc
import mudata as mu
import sys
import pandas as pd
from anndata import AnnData

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "output_key": "pca",
    "num_components": 25,
    "layer": None,
    "obsm_output": "X_pca",
    "var_input": "filter_with_hvg",
    "varm_output": "varm_output",
    "uns_output": "pca_variance",
    "overwrite": True,
}

meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading %s, modality %s", par["input"], par["modality"])
data = mu.read_h5ad(par["input"], mod=par["modality"])

logger.info("Computing PCA components for modality '%s'", par["modality"])
if par["layer"] and par["layer"] not in data.layers:
    raise ValueError(f"{par['layer']} was not found in modality {par['modality']}.")

chunked, chunk_size = par["chunked"], par["chunk_size"]
if chunked:
    if not chunk_size:
        raise ValueError(
            "Requested to perform an incremental PCA "
            "('chunked'), but the chunk size is not set."
        )
    if chunk_size < par["num_components"]:
        raise ValueError(
            f"The requested chunk size ({chunk_size}) must not be smaller "
            f"than the number of components ({par['num_components']})"
        )

layer = data.X if not par["layer"] else data.layers[par["layer"]]
adata_input_layer = AnnData(layer, var=pd.DataFrame([], index=data.var.index))

mask_var = None
if par["var_input"]:
    if par["var_input"] not in data.var.columns:
        raise ValueError(
            f"Requested to use .var column {par['var_input']} "
            "as a selection of genes to run the PCA on, "
            f"but the column is not available for modality {par['modality']}"
        )
    mask_var = data.var[par["var_input"]]

# run pca
sc.tl.pca(
    adata_input_layer,
    n_comps=par["num_components"],
    copy=False,  # A copy was already created
    return_info=True,
    mask_var=mask_var,
    chunked=chunked,
    chunk_size=chunk_size,
    random_state=par["seed"],
)

# store output in specific objects

check_exist_dict = {
    "obsm_output": ("obs"),
    "varm_output": ("varm"),
    "uns_output": ("uns"),
}
for parameter_name, field in check_exist_dict.items():
    if par[parameter_name] in getattr(data, field):
        if not par["overwrite"]:
            raise ValueError(
                f"Requested to create field {par[parameter_name]} in .{field} "
                f"for modality {par['modality']}, but field already exists."
            )
        del getattr(data, field)[par[parameter_name]]

data.obsm[par["obsm_output"]] = adata_input_layer.obsm["X_pca"]
data.varm[par["varm_output"]] = adata_input_layer.varm["PCs"]
data.uns[par["uns_output"]] = {
    "variance": adata_input_layer.uns["pca"]["variance"],
    "variance_ratio": adata_input_layer.uns["pca"]["variance_ratio"],
}


logger.info(
    "Writing to %s with compression %s.", par["output"], par["output_compression"]
)
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], data, par["output_compression"]
)

logger.info("Finished")
