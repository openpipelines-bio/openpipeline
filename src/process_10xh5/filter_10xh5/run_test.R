library(testthat)
requireNamespace("hdf5r", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
library(anndata)

cat("Run command\n")
input <- paste0(meta$resources_dir, "/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5")
output <- "out.h5"

system(paste0(
  "./", meta$functionality_name, " ",
  "--input '", input, "' ",
  "--output ", output, " ",
  "--min_library_size 1000 ",
  "--min_cells_per_gene 300 ",
  "--verbose"
))

sc <- reticulate::import("scanpy")
orig_ad <- sc$read_10x_h5(input)

out_ad <- sc$read_10x_h5(output)

# check if matrix is not empty
expect_gte(nrow(out_ad), 100)
expect_gte(ncol(out_ad), 100)

# check if filtering has been performed
expect_lte(nrow(out_ad), nrow(orig_ad))
expect_lte(ncol(out_ad), ncol(orig_ad))

cat(">> All tests passed successfully!\n")
