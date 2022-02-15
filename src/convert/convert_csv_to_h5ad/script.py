### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.csv",
    "output": "output.h5ad",
    "delimiter": ",",
    "use_column_name": "true",
    "compression": "gzip",
}
### VIASH END

import scanpy as sc
import scipy

print("Converting " + par["input"] + " to " + par["output"])

data = sc.read_csv(
    par["input"], delimiter=par["delimiter"], first_column_names=par["use_column_names"]
)

data.var_names_make_unique()
data.X = scipy.sparse.csr_matrix(data.X)
data.raw = data

data.write(par["output"], compression=par["compression"])
