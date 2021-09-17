### VIASH START

par ={
    "input": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.h5ad",
    "output": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.norm.hvg.h5ad"
}
### VIASH END

import argparse
import anndata
import scanpy as sc

data = anndata.read_h5ad(par["input"])

sc.pp.highly_variable_genes(data, flavor=par["flavor"])

if len(par["excludedGenes"]) > 0:
    excludedGenes = list(map(str.strip, par["excludedGenes"].split(",")))
    print("Excluding genes: " + str(excludedGenes))
    
    data.var["highly_variable"] = (~data.var["highly_variable"].index.isin(excludedGenes) & data.var["highly_variable"].values)

data.write(par["output"], compression="gzip")