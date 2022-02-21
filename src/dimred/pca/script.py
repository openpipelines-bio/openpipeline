import anndata
import scanpy as sc

### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.h5ad",
    "output": "output.h5ad",
    "nPCs": 25,
}
### VIASH END

print("Reading", par["input"])
data = anndata.read_h5ad(par["input"])

print("Computing PCA coordinates")
sc.tl.pca(data, n_comps=par["nPCs"])
data.uns["pcaParameters"] = {"PCA: nPCs": par["nPCs"]}

print("Writing".par["output"])
data.write(par["output"], compression="lzf")
