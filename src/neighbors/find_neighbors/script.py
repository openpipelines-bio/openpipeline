import sys
import numpy as np
numpy_module = sys.modules['numpy']
numpy_module.float_ = np.float64
sys.modules['numpy'] = numpy_module

import mudata as mu
import scanpy as sc
import sys

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "metric": 'cosine',
    "num_neighbors": 15,
    "modality": "rna",
    "obsm_input": "X_pca",
    "uns_output": "neighbors",
    "obsp_distances": "distances",
    "obsp_connectivities": "connectivities",
    "seed": None
}
meta = {
    'resources_dir': "."
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

logger.info("Reading input mudata")
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Computing a neighborhood graph on modality %s", mod)
adata = mdata.mod[mod]
neighbors = sc.Neighbors(adata)
neighbors.compute_neighbors(
    n_neighbors=par["num_neighbors"],
    use_rep=par["obsm_input"],
    metric=par["metric"],
    random_state=par["seed"],
    method="umap",
)

adata.uns[par["uns_output"]] = {
    'connectivities_key': par["obsp_connectivities"],
    'distances_key': par["obsp_distances"],
    'params': {
        'n_neighbors': neighbors.n_neighbors,
        'method': "umap",
        'random_state': par["seed"],
        'metric': par["metric"],
        'use_rep': par["obsm_input"]
    }
}

adata.obsp[par["obsp_distances"]] = neighbors.distances
adata.obsp[par["obsp_connectivities"]] = neighbors.connectivities


logger.info("Writing to %s", par["output"])
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])
