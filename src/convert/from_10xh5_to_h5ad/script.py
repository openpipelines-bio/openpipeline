import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5",
    "output": "foo.h5ad",
    "gex_only": False
}
## VIASH END

print("Reading", par["input"])
data = sc.read_10x_h5(par["input"], gex_only=par["gex_only"])

print("Making unique")
data.var_names_make_unique()

print("Writing", par["output"])
data.write_h5ad(par["output"])
