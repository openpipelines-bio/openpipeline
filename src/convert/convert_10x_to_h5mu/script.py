import muon as mu

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
    "output": "output.h5mu",
}
## VIASH END

print("Converting", par["input"], "to", par["output"])
mu.read_10x_h5(par["input"]).write(par["output"])
