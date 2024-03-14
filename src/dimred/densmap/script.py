from umap import UMAP
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

neigh_key = par["uns_neighbors"]
temp_uns = { neigh_key: data.uns[neigh_key] }
pca_key = temp_uns[neigh_key]['params']['use_rep']
knn_indices_key = temp_uns[neigh_key]['knn_indices_key']
knn_distances_key = temp_uns[neigh_key]['knn_distances_key']


X_densmap = UMAP(
  min_dist=par["min_dist"],
  spread=par["spread"],
  n_components=par["num_components"],
  n_epochs=par["max_iter"],
  learning_rate=par["alpha"],
  repulsion_strength=par["gamma"],
  negative_sample_rate=par["negative_sample_rate"],
  init=par["init_pos"],
  metric=data.uns["neighbors"]["metric"],
  metric_kwds=data.uns["neighbors"].get("metric_kwds", {}),
  densmap=True,
  dens_lambda=par["lambda"],
  dens_frac=par["fraction"],
  dens_var_shift=par["var_shift"],
  precomputed_knn=(
    data.obsm[knn_indices_key],
    data.obsm[knn_distances_key]
  )
).fit_transform(data.obsm[pca_key])

data.obsm[par['obsm_output']] = X_densmap

logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])

logger.info("Finished")