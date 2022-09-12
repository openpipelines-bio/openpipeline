import muon as mu
import scanpy as sc
import logging
from sys import stdout

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_ums.h5mu",
    "output": "output.h5mu",
    "metric": 'cosine',
    "num_neighbors": 15,
    "modality": ["rna"],
    "obsm_input": "X_pca",
    "obsp_output": "distances",
    "uns_output": "connectivities"
}
meta = {"functionality_name": "lognorm"}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])

for mod in par["modality"]:
    logger.info("Computing a neighborhood graph on modality %s", mod)
    adata = mdata.mod[mod]
    neighbors = sc.Neighbors(adata)
    neighbors.compute_neighbors(
        n_neighbors=par["num_neighbors"], 
        # use_rep=par["obsm_input"],
        metric=par["metric"]
    )
    adata.uns[par["uns_output"]] = neighbors.connectivities
    adata.obsp[par["obsp_output"]] = neighbors.distances

logger.info("Writing to %s", par["output"])
mdata.write(filename=par["output"])
