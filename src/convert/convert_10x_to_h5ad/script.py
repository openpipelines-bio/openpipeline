### VIASH START
par = {
	"input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
	"output": "test.h5ad",
	"compression": "gzip"
}
### VIASH END

import anndata
import scanpy as sc

print("Converting " + par["input"] + " to " + par["output"])
input_data = sc.read_10x_h5(par["input"], gex_only = False)

print("Setting Gene Expression as X")
adata = input_data[:, input_data.var["feature_types"] == "Gene Expression"] #.copy()

print("Making var_names unique")
try: 
	adata.var_names_make_unique()
except:
	pass

print("Setting raw")
adata.raw = adata

# adding more data to different obsm slots
def add_data_to_obsm(feature_type, obsm_field):
  more_data = input_data[:, input_data.var["feature_types"] == feature_type]
  
  if more_data.n_vars > 0:
    print("Saving " + feature_type + " data to adata.obsm[" + obsm_field + "]")
    adata.obsm[obsm_field] = more_data.to_df()

add_data_to_obsm("Antibody Capture", "counts_antibody")
add_data_to_obsm("CRISPR Guide Capture", "counts_crispr")
add_data_to_obsm("Custom", "counts_custom")

print("Writing adata file to " + par["output"])
adata.write(par["output"], compression = par["compression"])
