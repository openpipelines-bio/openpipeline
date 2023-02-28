import scanpy as sc
import mudata as mu
import numpy as np
import logging
from sys import stdout
import re

## VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
  'modality': 'rna',
  'output': 'output.h5mu',
  'var_name_filter': 'filter_with_hvg',
  'do_subset': False,
  'flavor': 'seurat',
  'n_top_genes': 123,
  'min_mean': 0.0125,
  'max_mean': 3.0,
  'min_disp': 0.5,
  'span': 0.3,
  'n_bins': 20,
  'varm_name': 'hvg',
  'layer': 'log_normalized'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

mdata = mu.read_h5mu(par["input"])
mdata.var_names_make_unique()

mod = par['modality']
logger.info(f"Processing modality '%s'", mod)
data = mdata.mod[mod]

# Workaround for issue
# https://github.com/scverse/scanpy/issues/2239
# https://github.com/scverse/scanpy/issues/2181
if par['flavor'] != "seurat_v3":
    # This component requires log normalized data when flavor is not seurat_v3
    # We assume that the data is correctly normalized but scanpy will look at
    # .uns to check the transformations performed on the data.
    # To prevent scanpy from automatically tranforming the counts when they are
    # already transformed, we set the appropriate values to .uns.
    if 'log1p' not in data.uns:
        logger.warning("When flavor is not set to 'seurat_v3', "
                       "the input data for this component must be log-transformed. "
                       "However, the 'log1p' dictionairy in .uns has not been set. "
                       "This is fine if you did not log transform your data with scanpy."
                       "Otherwise, please check if you are providing log transformed "
                       "data using --layer.")
        data.uns['log1p'] = {'base': None}
    elif 'log1p' in data.uns and 'base' not in data.uns['log1p']:
        data.uns['log1p']['base'] = None

logger.info("\tUnfiltered data: %s", data)

logger.info("\tComputing hvg")
# construct arguments
hvg_args = {
    'adata': data,
    'n_top_genes': par["n_top_genes"],
    'min_mean': par["min_mean"],
    'max_mean': par["max_mean"],
    'min_disp': par["min_disp"],
    'span': par["span"],
    'n_bins': par["n_bins"],
    'flavor': par["flavor"],
    'subset': False,
    'inplace': False,
    'layer': par['layer'],
}

optional_parameters = {
    "max_disp": "max_disp",
    "obs_batch_key": "batch_key",
    "n_top_genes": "n_top_genes"
}
# only add parameter if it's passed
for par_name, dest_name in optional_parameters.items():
    if par.get(par_name):
        hvg_args[dest_name] = par[par_name]

# scanpy does not do this check, although it is stated in the documentation
if par['flavor'] == "seurat_v3" and not par['n_top_genes']:
    raise ValueError("When flavor is set to 'seurat_v3', you are required to set 'n_top_genes'.")

# call function
try:
    out = sc.pp.highly_variable_genes(**hvg_args)
    out.index = data.var.index
except ValueError as err:
    if str(err) == "cannot specify integer `bins` when input data contains infinity":
        err.args = ("Cannot specify integer `bins` when input data contains infinity. "
                    "Perhaps input data has not been log normalized?",)
    if re.search("Bin edges must be unique:", str(err)):
        raise RuntimeError("Scanpy failed to calculate hvg. The error "
                           "returned by scanpy (see above) could be the "
                           "result from trying to use this component on unfiltered data.") from err
    raise err

logger.info("\tStoring output into .var")
if par.get("var_name_filter", None) is not None:
    data.var[par["var_name_filter"]] = out["highly_variable"]

if par.get("varm_name", None) is not None and 'mean_bin' in out:
    # drop mean_bin as mudata/anndata doesn't support tuples
    data.varm[par["varm_name"]] = out.drop("mean_bin", axis=1)

if par["do_subset"]:
    keep_feats = np.ravel(data.var[par["var_name_filter"]])
    mdata.mod[mod] = data[:,keep_feats]

logger.info("Writing h5mu to file")
mdata.write_h5mu(par["output"], compression="gzip")
