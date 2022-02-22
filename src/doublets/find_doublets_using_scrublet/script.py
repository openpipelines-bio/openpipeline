import scrublet as scr
import scanpy as sc
import muon as mu

### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "expected_doublet_rate": 0.05,
    "min_counts": 2,
    "min_cells": 3,
    "min_gene_variablity_percent": 85,
    "prinipal_components_amount": 30,
    "distance_metric": "euclidean",
    "col_name_doublet_score": "scrublet_doublet_score",
    "col_name_predicted_doublets": "scrublet_predicted_doublets",
    "outputFormat": "h5mu",
}
### VIASH END

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
