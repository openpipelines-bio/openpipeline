library(MuDataSeurat)

## VIASH START
par <- list(
  input = "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
  output = "output.rds"
)
## VIASH END

cat("Reading input file\n")
obj <- ReadH5MU(par$input)

cat("Writing output file\n")
saveRDS(obj, file = par$output, compress = TRUE)