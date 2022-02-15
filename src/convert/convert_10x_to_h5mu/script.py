### VIASH START
par = {
    "input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
    "output": "test.h5mu",
}
### VIASH END

import muon as mu

print("Converting", par["input"], "to", par["output"])
mu.read_10x_h5(par["input"]).write(par["output"])
