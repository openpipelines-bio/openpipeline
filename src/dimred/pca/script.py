### VIASH START

par ={
    "input": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.h5ad",
    "output": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.pca.h5ad",
    "nPCs": 20
}
### VIASH END

import argparse
import anndata
import scanpy as sc

data = anndata.read_h5ad(par["input"])

        
sc.tl.pca(data, 
          n_comps = par["nPCs"])

data.uns["pcaParameters"] = {
    "PCA: nPCs": par["nPCs"]
}


data.write(par["output"], compression = "lzf")