### VIASH START

par = {
    "input": "./test/",
    "output": "./test/",
    
    "expectedDoubletRate": 0.06,
    "minCounts": 2,
    "minCells": 3, 
    "minGeneVariabilityPercentile": 85,
    "nPrincipalComponents": 30,
    "distanceMetric": "euclidean",
    "colNameDoubletScore": "ScrubletDoubletScore",
    "colNamePredictedDoublets": "ScrubletPredictedDoublets"
}


### VIASH END

import scrublet as scr
import scanpy as sc

data = sc.read_h5ad(par["input"])

scrub = scr.Scrublet(data.X)

doubletScores, predictedDoublets = scrub.scrub_doublets(min_counts=par["minCounts"], 
                                                          min_cells=par["minCells"], 
                                                          min_gene_variability_pctl=par["minGeneVariabilityPercentile"], 
                                                          n_prin_comps=par["nPrincipalComponents"],
                                                          distance_metric=par["distanceMetric"])


data.obs[par["colNameDoubletScore"]] = doubletScores
data.obs[par["colNamePredictedDoublets"]] = predictedDoublets

if (par["outputFormat"] == "h5ad"):
    data.obs[par["colNameDoubletScore"]] = doubletScores
    data.obs[par["colNamePredictedDoublets"]] = predictedDoublets

    data.write_h5ad(par["output"], compression = "gzip")

elif (par["outputFormat"] == "csv"):
     data.obs[[par["colNameDoubletScore"], par["colNamePredictedDoublets"]]].to_csv(par["output"])

else: 
     raise ValueError("An unrecognized output format was specified.")