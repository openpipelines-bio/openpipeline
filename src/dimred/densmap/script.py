from umap.umap_ import find_ab_params, simplicial_set_embedding
from sklearn.utils import check_random_state
import mudata as mu
import sys
import anndata as ad

## VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu',
  'modality': 'rna',
  'output': 'output.h5mu',
  'obsm_output': 'X_densmap',
  'lambda': 2.0,
  'fraction': 0.3,
  'var_shift': 0.1
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

logger.info("Reading %s", par["input"])
mdata = mu.read_h5mu(par["input"])

logger.info("Computing densMAP for modality '%s'", par['modality'])
data = mdata.mod[par['modality']]

if par['uns_neighbors'] not in data.uns:
    raise ValueError(f"'{par['uns_neighbors']}' was not found in .mod['{par['modality']}'].uns.")

# create temporary AnnData
# ... because sc.tl.umap doesn't allow to choose
# the obsm output slot
# ... also we can see scanpy is a data format dependency hell
neigh_key = par["uns_neighbors"]
temp_uns = { neigh_key: data.uns[neigh_key] }
conn_key = temp_uns[neigh_key]['connectivities_key']
dist_key = temp_uns[neigh_key]['distances_key']
temp_obsp = {
  conn_key: data.obsp[conn_key],
  dist_key: data.obsp[dist_key],
}
pca_key = temp_uns[neigh_key]['params']['use_rep']
temp_obsm = {
  pca_key: data.obsm[pca_key]
}

default_epochs = 500 if data.obsp[conn_key].shape[0] <= 10000 else 200
n_epochs = default_epochs if par["max_iter"] is None else par["max_iter"]
a, b = find_ab_params(par["spread"], par["mim_dist"])

X_densmap = simplicial_set_embedding(
  data=data.obsm[pca_key],
  graph=data.obsp[conn_key],
  n_components=par["n_components"],
  initial_alpha=par["alpha"],
  a=a,
  b=b,
  gamma=par["gamma"],
  negative_sample_rate=par["negative_sample_rate"],
  n_epochs=par["max_iter"],
  init=par["init_pos"],
  random_state = check_random_state(None),
  metric=data.uns[neigh_key].get("metric", "euclidean"),
  metric_kwds=data.uns[neigh_key].get("metric_kwds", {}),
  densmap=True,
  densmap_kwds={
    "lambda": par["lambda"],
    "fraction": par["fraction"],
    "var_shift": par["var_shift"]
  }
)