### VIASH START

par ={
    "input": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.h5ad",
    "output": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.h5ad",
    "metric": "cosine",
    "nNeighbors": 30
}
### VIASH END

import anndata

data = anndata.read_h5ad(par["input"])

sc.pp.neighbors(data, 
    n_neighbors = par["nNeighbors"], 
    metric = par["metric"])

data.uns["nearestNeighbourParameters"] = {
    "Nearest neighbours: metric": par["metric"],
    "Nearest neighbours: number of neighbors": par["nNeighbors"]
}

data.write(par["output"], compression = par["compression"])
