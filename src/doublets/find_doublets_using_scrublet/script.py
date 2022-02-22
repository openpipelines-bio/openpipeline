### VIASH START

par = {
    "input": "./test/",
    "output": "./test/",
    "expected_doublet_rate": 0.06,
    "min_counts": 2,
    "min_cells": 3,
    "min_gene_variablity_percent": 85,
    "prinipal_components_amount": 30,
    "distance_metric": "euclidean",
    "col_name_doublet_score": "scrublet_doublet_score",
    "col_name_predicted_doublets": "scrublet_predicted_doublets",
}


### VIASH END

import scrublet as scr
import scanpy as sc
import muon as mu

# data = sc.read_h5ad(par["input"])

print("Reading", par["input"])
mdata = mu.read_h5mu(par["input"])

scrub = scr.Scrublet(mdata.mod["rna"].X)

doubletScores, predictedDoublets = scrub.scrub_doublets(
    min_counts=par["min_counts"],
    min_cells=par["min_cells"],
    min_gene_variability_pctl=par["min_gene_variablity_percent"],
    n_prin_comps=par["prinipal_components_amount"],
    distance_metric=par["distance_metric"],
)

mdata.obs[par["col_name_doublet_score"]] = doubletScores
mdata.obs[par["col_name_predicted_doublets"]] = predictedDoublets

if par["outputFormat"] == "h5mu":
    mdata.obs[par["col_name_doublet_score"]] = doubletScores
    mdata.obs[par["col_name_predicted_doublets"]] = predictedDoublets

    print("Writing", par["output"])
    mdata.write_h5mu(filename=par["output"])
else:
    print(par["outputFormat"], " format is not supported. No file will be written.")
