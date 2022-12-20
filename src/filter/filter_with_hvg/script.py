import scanpy as sc
import mudata as mu
import numpy as np
import logging
from sys import stdout

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
if 'log1p' in data.uns and 'base' not in data.uns['log1p']:
    data.uns['log1p']['base'] = None
#sc.pp.log1p(data)

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
    'layer': par['layer']
}

# only add parameter if it's passed
if par.get("max_disp", None) is not None:
    hvg_args["max_disp"] = par["max_disp"]
if par.get("obs_batch_key", None) is not None:
    hvg_args["batch_key"] = par["obs_batch_key"]

# call function
try:
    out = sc.pp.highly_variable_genes(**hvg_args)
    out.index = data.var.index
except ValueError as err:
    if str(err) == "cannot specify integer `bins` when input data contains infinity":
        err.args = ("Cannot specify integer `bins` when input data contains infinity. Perhaps input data has not been log normalized?",)
    raise err

logger.info("\tStoring output into .var")
if par.get("var_name_filter", None) is not None:
    data.var[par["var_name_filter"]] = out["highly_variable"]

if par.get("varm_name", None) is not None:
    # drop mean_bin as mudata/anndata doesn't support tuples
    data.varm[par["varm_name"]] = out.drop("mean_bin", axis=1)

if par["do_subset"]:
    keep_feats = np.ravel(data.var[par["var_name_filter"]])
    mdata.mod[mod] = data[:,keep_feats]

logger.info("Writing h5mu to file")
mdata.write_h5mu(par["output"], compression="gzip")
