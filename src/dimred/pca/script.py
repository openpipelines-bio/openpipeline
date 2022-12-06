import scanpy as sc
import mudata as mu
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "output_key": "pca",
    "num_components": 25
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading %s.", par["input"])
mdata = mu.read_h5mu(par["input"])

logger.info("Computing PCA components for modality '%s'", par['modality'])
data = mdata.mod[par['modality']]
if par['layer'] and par['layer'] not in data.layers:
    raise ValueError(f"{par['layer']} was not found in modality {par['modality']}.")
layer = data.X if not par['layer'] else data.layers[par['layer']]
# run pca
X_pca, loadings, variance, variance_ratio = sc.tl.pca(
    layer, 
    n_comps=par["num_components"], 
    return_info=True
)

# store output in specific objects
data.obsm[par["obsm_output"]] = X_pca
data.varm[par["varm_output"]] = loadings.T
data.uns[par["uns_output"]] = { "variance": variance, "variance_ratio": variance_ratio }

logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"])

logger.info("Finished")