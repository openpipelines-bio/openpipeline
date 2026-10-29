import scanpy as sc
import mudata as mu
import anndata as ad
import numpy as np
import pandas as pd
import sys
import re

ad.settings.allow_write_nullable_strings = True


## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "var_name_filter": "filter_with_hvg",
    "do_subset": False,
    "flavor": "seurat_v3",
    "n_top_features": 20,
    "min_mean": 0.0125,
    "max_mean": 3.0,
    "min_disp": 0.5,
    "span": 0.3,
    "n_bins": 20,
    "var_input": None,
    "features_to_exclude": ["ENSG00000237613"],
    "output_compression": "gzip",
    "varm_name": "hvg",
    "obs_batch_key": "batch",
    "layer": "log_transformed",
}

meta = {"resources_dir": "src/utils/"}

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
column_index = rna_in.var.columns.get_indexer(["common_vars"])
rna_in.var.iloc[:10000, column_index] = True
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
input_anndata = ad.AnnData(X=input_layer.copy(), obs=obs, var=var)
if "log1p" in data.uns:
    input_anndata.uns["log1p"] = data.uns["log1p"]

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

# Enable calculating the HVG only on a subset of vars
# e.g for cell type annotation, only calculate HVG on variables that are common between query and reference
if par["var_input"]:
    input_anndata.var[par["var_input"]] = data.var[par["var_input"]]
    input_anndata = subset_vars(input_anndata, par["var_input"])

# Exclude user-specified features from HVG calculation
excluded_features_mask = None
if par.get("features_to_exclude"):
    features_to_exclude = set(par["features_to_exclude"])
    logger.info(
        "\tExcluding %d specified features from HVG calculation",
        len(features_to_exclude),
    )
    excluded_features_mask = input_anndata.var_names.isin(features_to_exclude)
    n_excluded = excluded_features_mask.sum()
    n_not_found = len(features_to_exclude) - n_excluded
    if n_not_found > 0:
        not_found = features_to_exclude - set(
            input_anndata.var_names[excluded_features_mask]
        )
        logger.warning(
            "\t%d features to exclude were not found in the data: %s",
            n_not_found,
            list(not_found)[:10],
        )
    logger.info("\tExcluding %d features from HVG calculation", n_excluded)
    if n_excluded == input_anndata.n_vars:
        raise ValueError(
            f"All features ({n_excluded}) are in the exclusion list. "
            "Please check your --features_to_exclude list."
        )
    # Subset to non-excluded features for HVG calculation using subset_vars
    input_anndata = subset_vars(input_anndata, ~excluded_features_mask)
    logger.info("\t%d features remaining for HVG calculation", input_anndata.n_vars)

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


def align_output_to_var(df_with_missing_elements, target_var):
    # Make sure string columns become a nullable dtype
    df_with_missing_elements = df_with_missing_elements.convert_dtypes(
        infer_objects=True,
        convert_string=True,
        convert_integer=False,
        convert_boolean=False,
        convert_floating=False,
    )
    # The reindex below matches rows by feature name. If the index of the output
    # does not hold the feature names, every row would silently become NA
    # (and all features would be marked as not highly variable).
    unknown_features = ~df_with_missing_elements.index.isin(target_var.index)
    assert not unknown_features.any(), (
        f"'highly_variable_genes' output contains {unknown_features.sum()} "
        "features that are not present in the input, expected the index to "
        "contain the feature names."
    )

    fill_vals = {
        "means": np.nan,
        "gene_name": pd.NA,
        "mean_bin": np.nan,
        "highly_variable": False,
        "dispersions": np.nan,
        "dispersions_norm": np.nan,
        "variances": np.nan,
        "variances_norm": np.nan,
        "highly_variable_rank": np.nan,
        "highly_variable_nbatches": np.nan,
        "highly_variable_intersection": False,
    }
    unexpected_columns = df_with_missing_elements.columns.difference(fill_vals.keys())
    if not unexpected_columns.empty:
        raise RuntimeError(
            f"'highly_variable_genes' output contains unexpected columns: {''.join(unexpected_columns.to_list())}"
        )

    # Reindex each column separately with its own fill value. Filling while
    # reindexing keeps the dtype (e.g. bool), while reindexing first and
    # calling fillna afterwards would turn bool columns into object columns.
    return pd.DataFrame(
        {
            column: values.reindex(target_var.index, fill_value=fill_vals[column])
            for column, values in df_with_missing_elements.items()
        },
        index=target_var.index,
    )


# call function
try:
    out = sc.pp.highly_variable_genes(**hvg_args)
    out = align_output_to_var(out, data.var)
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

logger.info("\tStoring output into .var")
if par.get("var_name_filter", None) is not None:
    data.var[par["var_name_filter"]] = out["highly_variable"]

if par.get("varm_name", None) is not None:
    # drop mean_bin if present as mudata/anndata doesn't support tuples
    if "mean_bin" in out.columns:
        out = out.drop(columns=["mean_bin"])
    data.varm[par["varm_name"]] = out

if par.get("uns_name", None) is not None:
    data.uns[par["uns_name"]] = {"flavor": par["flavor"]}

logger.info("Writing h5mu to file")
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], modality_name, data, par["output_compression"]
)
