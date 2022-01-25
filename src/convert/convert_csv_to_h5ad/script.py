### VIASH START
par = {
	"input": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.csv",
	"output": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad"
}
### VIASH END

import scanpy as sc
import scipy

print("Converting " + par["input"] + " to " + par["output"])

data = sc.read_csv(par["input"], 
	delimiter = par["delimiter"],
	first_column_names = par["useColumnNames"])   

data.var_names_make_unique()
data.X = scipy.sparse.csr_matrix(data.X)
data.raw = data

data.write(par["output"], compression = par["compression"])