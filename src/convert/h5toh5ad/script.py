### VIASH START

par = {
	"input": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
	"output": "./test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad"
}
### VIASH END


import argparse
import anndata
import scanpy as sc


print("Converting " + par["input"] + " to " + par["output"])

data = sc.read_10x_h5(par["input"], gex_only = False)       

gexData = data[:, data.var["feature_types"] == "Gene Expression"]

try: 
	gexData.var_names_make_unique()
except:
	pass
gexData.raw = gexData

def addDataToObsm(featureType, obsmField):
	nData = data[:, data.var["feature_types"] == featureType].to_df()
	gexData.obsm[obsmField] = nData

addDataToObsm("Antibody Capture", "counts_antibody")
addDataToObsm("CRISPR Guide Capture", "counts_crispr")
addDataToObsm("Custom", "counts_custom")

gexData.write(par["output"], compression = par["compression"])
