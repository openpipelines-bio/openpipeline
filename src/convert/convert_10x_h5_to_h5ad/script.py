

### VIASH START

par = {
	"input": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
	"output": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad"
}
### VIASH END

import scanpy as sc

data = sc.read_10x_h5(par["input"], gex_only = False)       

d = data[:, data.var["feature_types"] == par["modality"]]

try: 
	d.var_names_make_unique()
except:
	pass

d.raw = d

d.write(par["output"], compression = par["compression"])
