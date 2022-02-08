### VIASH START
par = {
	"input": "resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5",
	"output": "test.h5mu",
	"compression": "gzip"
}
### VIASH END

import anndata
import muon as mu

mu.read_10x_h5(par["input"]).write(par["output"], compression = par["compression"])
