import scanpy as sc
import mudata as mu
import anndata as ad
import pandas as pd
import sys
import re

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "var_name_filter": "filter_with_hvg",
    "do_subset": False,
    "flavor": "seurat",
    "n_top_features": None,
    "min_mean": 0.0125,
    "max_mean": 3.0,
    "min_disp": 0.5,
    "span": 0.3,
    "n_bins": 20,
    "varm_name": "hvg",
    "obs_batch_key": "batch",
    "layer": "log_transformed",
    "var_input": "common_vars",
}

meta = {"resources_dir": "src/utils"}

mu_in = mu.read_h5mu(
    "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)
rna_in = mu_in.mod["rna"]
assert "filter_with_hvg" not in rna_in.var.columns
log_transformed = sc.pp.log1p(rna_in, copy=True)
rna_in.layers["log_transformed"] = log_transformed.X
rna_in.uns["log1p"] = log_transformed.uns["log1p"]
temp_h5mu = "lognormed.h5mu"
rna_in.obs["batch"] = "A"
column_index = rna_in.obs.columns.get_indexer(["batch"])
rna_in.obs.iloc[slice(rna_in.n_obs // 2, None), column_index] = "B"
rna_in.var["common_vars"] = False
rna_in.var["common_vars"].iloc[:10000] = True
mu_in.write_h5mu(temp_h5mu)
par["input"] = temp_h5mu
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from subset_vars import subset_vars
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()
modality_name = par["modality"]
data = mu.read_h5ad(par["input"], mod=modality_name)
assert data.var_names.is_unique, "Expected var_names of input modality to be unique."

logger.info("Processing modality '%s'", modality_name)

if par["layer"] and par["layer"] not in data.layers:
    raise ValueError(
        f"Layer '{par['layer']}' not found in layers for modality '{modality_name}'. "
        f"Found layers are: {','.join(data.layers)}"
    )

# input layer argument does not work when batch_key is specified because
# it still uses .X to filter out genes with 0 counts, even if .X might not exist.
# So create a custom anndata as input that always uses .X
input_layer = data.X if not par["layer"] else data.layers[par["layer"]]
obs = pd.DataFrame(index=data.obs_names.copy())
var = pd.DataFrame(index=data.var_names.copy())
if par["obs_batch_key"]:
    obs = data.obs.loc[:, par["obs_batch_key"]].to_frame()
if par["var_input"]:
    var = data.var.loc[:, par["var_input"]].to_frame()
input_anndata = ad.AnnData(X=input_layer.copy(), obs=obs, var=var)
if "log1p" in data.uns:
    input_anndata.uns["log1p"] = data.uns["log1p"]

# Workaround for issue
# https://github.com/scverse/scanpy/issues/2239
# https://github.com/scverse/scanpy/issues/2181
if par["flavor"] != "seurat_v3":
    # This component requires log normalized data when flavor is not seurat_v3
    # We assume that the data is correctly normalized but scanpy will look at
    # .uns to check the transformations performed on the data.
    # To prevent scanpy from automatically tranforming the counts when they are
    # already transformed, we set the appropriate values to .uns.
    if "log1p" not in input_anndata.uns:
        logger.warning(
            "When flavor is not set to 'seurat_v3', "
            "the input data for this component must be log-transformed. "
            "However, the 'log1p' dictionairy in .uns has not been set. "
            "This is fine if you did not log transform your data with scanpy."
            "Otherwise, please check if you are providing log transformed "
            "data using --layer."
        )
        input_anndata.uns["log1p"] = {"base": None}
    elif "log1p" in input_anndata.uns and "base" not in input_anndata.uns["log1p"]:
        input_anndata.uns["log1p"]["base"] = None

# Enable calculating the HVG only on a subset of vars
# e.g for cell type annotation, only calculat HVG on vars that are common between query and reference
if par["var_input"]:
    input_anndata = subset_vars(input_anndata, par["var_input"])

logger.info("\tUnfiltered data: %s", data)

logger.info("\tComputing hvg")
# construct arguments
hvg_args = {
    "adata": input_anndata,
    "n_top_genes": par["n_top_features"],
    "min_mean": par["min_mean"],
    "max_mean": par["max_mean"],
    "min_disp": par["min_disp"],
    "span": par["span"],
    "n_bins": par["n_bins"],
    "flavor": par["flavor"],
    "subset": False,
    "inplace": False,
    "layer": None,  # Always uses .X because the input layer was already handled
}

optional_parameters = {
    "max_disp": "max_disp",
    "obs_batch_key": "batch_key",
    "n_top_genes": "n_top_features",
}
# only add parameter if it's passed
for par_name, dest_name in optional_parameters.items():
    if par.get(par_name):
        hvg_args[dest_name] = par[par_name]

# scanpy does not do this check, although it is stated in the documentation
if par["flavor"] == "seurat_v3" and not par["n_top_features"]:
    raise ValueError(
        "When flavor is set to 'seurat_v3', you are required to set 'n_top_features'."
    )

# call function
try:
    out = sc.pp.highly_variable_genes(**hvg_args)
    if par["obs_batch_key"] is not None:
        out = out.reindex(index=data.var.index, method=None)
        assert (
            out.index == data.var.index
        ).all(), "Expected output index values to be equivalent to the input index"
except ValueError as err:
    if str(err) == "cannot specify integer `bins` when input data contains infinity":
        err.args = (
            "Cannot specify integer `bins` when input data contains infinity. "
            "Perhaps input data has not been log normalized?",
        )
    if re.search("Bin edges must be unique:", str(err)):
        raise RuntimeError(
            "Scanpy failed to calculate hvg. The error "
            "returned by scanpy (see above) could be the "
            "result from trying to use this component on unfiltered data."
        ) from err
    raise err

out.index = data.var.index
logger.info("\tStoring output into .var")
if par.get("var_name_filter", None) is not None:
    data.var[par["var_name_filter"]] = out["highly_variable"]

if par.get("varm_name", None) is not None and "mean_bin" in out:
    # drop mean_bin as mudata/anndata doesn't support tuples
    data.varm[par["varm_name"]] = out.drop("mean_bin", axis=1)

logger.info("Writing h5mu to file")
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], modality_name, data, par["output_compression"]
)
