import scrublet as scr
import mudata as mu
import numpy as np
import logging
from sys import stdout

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

### VIASH START
par = {
    # "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "input": "output/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu.h5mu",
    "modality": [ "rna" ],
    "output": "output.h5mu",
    "obs_name_filter": "filter_with_scrublet",
    # "expected_doublet_rate": 0.05,
    "min_counts": 2,
    "min_cells": 3,
    "min_gene_variablity_percent": 85,
    "num_pca_components": 30,
    "distance_metric": "euclidean",
    "obs_name_doublet_score": "scrublet_doublet_score",
    "obs_name_predicted_doublets": "scrublet_predicted_doublets",
    "do_subset": True
}
meta = {
    'functionality_name': 'scrublet'
}
### VIASH END

logger.info("Reading %s.", par['input'])
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Processing modality '%s'.", mod)
data = mdata.mod[mod]

logger.info("\tRunning scrublet")
scrub = scr.Scrublet(data.X)

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=par["min_counts"],
    min_cells=par["min_cells"],
    min_gene_variability_pctl=par["min_gene_variablity_percent"],
    n_prin_comps=par["num_pca_components"],
    distance_metric=par["distance_metric"],
    use_approx_neighbors=False
)
keep_cells = np.invert(predicted_doublets)

logger.info("\tStoring output into .obs")
if par["obs_name_doublet_score"] is not None:
    data.obs[par["obs_name_doublet_score"]] = doublet_scores
if par["obs_name_filter"] is not None:
    data.obs[par["obs_name_filter"]] = keep_cells

if par["do_subset"]:
    mdata.mod[mod] = data[keep_cells, :]

logger.info("Writing h5mu to %s", par["output"])
mdata.write_h5mu(par["output"], compression="gzip")
