import scanpy as sc

### VIASH START
par = {
	"input": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
	"output": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",
	"modality": "Gene Expression",
	"compression": "gzip",
}
### VIASH END

data = sc.read_10x_h5(par["input"], gex_only=False)       


out = data[:, data.var["feature_types"] == par["modality"]]

try: 
	out.var_names_make_unique()
except:
	pass

out.raw = out

out.write_h5ad(par["output"], compression=par["compression"])
