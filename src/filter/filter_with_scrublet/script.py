import scrublet as scr
import mudata as mu
import numpy as np
import sys
import pandas as pd

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    # "input": "output/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "output_compression": "gzip",
    "obs_name_filter": "filter_with_scrublet",
    # "expected_doublet_rate": 0.05,
    "min_counts": 2,
    "min_cells": 3,
    "min_gene_variablity_percent": 85,
    "num_pca_components": 30,
    "distance_metric": "euclidean",
    "obs_name_doublet_score": "scrublet_doublet_score",
    "obs_name_predicted_doublets": "scrublet_predicted_doublets",
    "do_subset": True,
    "layer": None,
}
meta = {
    'name': 'scrublet',
    'resources_dir': '.',
}
### VIASH END

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

logger.info("Reading %s.", par['input'])
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Processing modality '%s'.", mod)
data = mdata.mod[mod]

logger.info("Using layer '%s'.", "X" if not par["layer"] else par["layer"])
input_layer = data.X if not par["layer"] else data.layers[par["layer"]]

logger.info("\tRunning scrublet")
scrub = scr.Scrublet(input_layer)

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=par["min_counts"],
    min_cells=par["min_cells"],
    min_gene_variability_pctl=par["min_gene_variablity_percent"],
    n_prin_comps=par["num_pca_components"],
    distance_metric=par["distance_metric"],
    use_approx_neighbors=False
)

try:
    keep_cells = np.invert(predicted_doublets)
except TypeError:
    if par['allow_automatic_threshold_detection_fail']:
        # Scrublet might not throw an error and return None if it fails to detect doublets...
        logger.info("\tScrublet could not automatically detect the doublet score threshold. Setting output columns to NA.")
        keep_cells = np.nan
        doublet_scores = np.nan
    else:
        raise RuntimeError("Scrublet could not automatically detect the doublet score threshold. "
                           "--allow_automatic_threshold_detection_fail can be used to ignore this failure "
                           "and set the corresponding output columns to NA.")

logger.info("\tStoring output into .obs")
if par["obs_name_doublet_score"] is not None:
    data.obs[par["obs_name_doublet_score"]] = doublet_scores
    data.obs[par["obs_name_doublet_score"]] = data.obs[par["obs_name_doublet_score"]].astype("float64")
if par["obs_name_filter"] is not None:
    data.obs[par["obs_name_filter"]] = keep_cells
    data.obs[par["obs_name_filter"]] = data.obs[par["obs_name_filter"]].astype(pd.BooleanDtype())

if par["do_subset"]:
    if pd.api.types.is_scalar(keep_cells) and pd.isna(keep_cells):
        logger.warning("Not subsetting beacuse doublets were not predicted")
    else: 
        mdata.mod[mod] = data[keep_cells, :]

logger.info("Writing h5mu to %s", par["output"])
mdata.write_h5mu(par["output"], compression=par["output_compression"])
