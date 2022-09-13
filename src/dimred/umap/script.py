import scanpy as sc
import muon as mu
import logging
from sys import stdout

## VIASH START
par = {
  'input': 'work/ca/588b7fcdfd953a534a59d794671451/pbmc_1k_protein_v3_mms.find_neighbors.output.h5mu',
  'modality': ['rna'],
  'output': 'output.h5mu',
  'output_key': 'umap',
  'min_dist': 0.5,
  'spread': 1.0,
  'num_components': 2,
  'max_iter': None,
  'alpha': 1.0,
  'gamma': 1.0,
  'negative_sample_rate': 5,
  'init_pos': 'spectral',
  'uns_neighbors': 'neighbors'
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading %s", par["input"])
mdata = mu.read_h5mu(par["input"])

for mod in par['modality']:
    logger.info("Computing UMAP for modality '%s'", mod)
    data = mdata.mod[mod]

    sc.tl.umap(
        data,
        min_dist=par["min_dist"],
        spread=par["spread"],
        n_components=par["num_components"],
        maxiter=par["max_iter"],
        alpha=par["alpha"],
        gamma=par["gamma"],
        negative_sample_rate=par["negative_sample_rate"],
        init_pos=par["init_pos"],
        neighbors_key=par["uns_neighbors"]
    )
    # note: should be able to set the neighbors key

logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"])

logger.info("Finished")