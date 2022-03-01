import muon as mu
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.nn.umap.h5mu",
    "output": "output.h5mu",
    'metric': 'cosine',
    'nNeighbors': 15,
}
## VIASH END

print("Reading", par["input"])
mdata = mu.read_h5mu(par["input"])

print("Compute a neighborhood graph of observations")
sc.pp.neighbors(mdata.mod["rna"], 
    n_neighbors = par["nNeighbors"], 
    metric = par["metric"])

mdata.uns["nearestNeighbourParameters"] = {
    "Nearest neighbours: metric": par["metric"],
    "Nearest neighbours: number of neighbors": par["nNeighbors"]
}

print("Writing", par["output"])
mdata.write_h5mu(filename=par["output"])
