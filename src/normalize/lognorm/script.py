### VIASH START

par = {
	"input": "./pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",
	"output": "./pbmc_1k_protein_v3_filtered_feature_bc_matrix.normalize.h5ad",

	"normalizedUMICount": "10000",
	"regressOutVariables": []
}

### VIASH END

import anndata
import scanpy as sc
from itertools import compress
import numpy as np
import multiprocessing


data = anndata.read_h5ad(par["input"])

print("Performing log normalization ... ")
sc.pp.normalize_total(data, target_sum = par["normalizedUMICount"])
sc.pp.log1p(data)

if any(map(len, par["regressOutVariables"])) > 0:
    selectNonEmpty = [len(i) > 0 for i in par["regressOutVariables"]]
    regressOutVariables = list(compress(par["regressOutVariables"], selectNonEmpty))
    
    sc.pp.regress_out(data, regressOutVariables, n_jobs=multiprocessing.cpu_count() - 1)

    data.uns["normalizationParameters"] = {
        "Normalization: method": "lognorm",
        "Normalization: normalizedUMICount": par["normalizedUMICount"],
        "Normalization: regressOutVariables": regressOutVariables
    }
else:
    data.uns["normalizationParameters"] = {
        "Normalization: method": "lognorm",
        "Normalization: normalizedUMICount": par["normalizedUMICount"]
    }
                
data.write(par["output"], compression = "lzf")       
    
