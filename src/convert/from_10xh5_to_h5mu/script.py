import muon as mu

## VIASH START
par = {
    "input": "5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5",
    "output": "5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5mu",
}
## VIASH END

print("Reading", par["input"])
mdata = mu.read_10x_h5(par["input"])

print("Making unique")
mdata.var_names_make_unique()

print("Writing", par["output"])
mdata.write_h5mu(par["output"])
