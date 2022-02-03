### VIASH START
par = {}
### VIASH END

import scanpy as sc
import scipy

data = sc.read_10x_mtx(par["input"])        

data.var_names_make_unique()

data = data[:, data.var["feature_types"] == par["modality"]]

data.X = scipy.sparse.csr_matrix(data.X)

data.raw = data

data.write(par["output"], compression = par["compression"])
    
