library(MuDataSeurat)

## VIASH START
par <- list(
  input = "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tx_processing.h5mu",
  output = "output.rds"
)
## VIASH END

obj <- ReadH5MU(par$input)

saveRDS(obj, file = par$output, compress = TRUE)