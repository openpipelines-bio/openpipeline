import muon as mu

## VIASH START
par = {
    "input": "5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix",
    "output": "5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5mu",
}
## VIASH END

mdata = mu.read_10x_mtx(par["input"])
mdata.write_h5mu(filename=par["output"])
