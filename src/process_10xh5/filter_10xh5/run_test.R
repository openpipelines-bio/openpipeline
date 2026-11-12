library(testthat)
requireNamespace("hdf5r", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
library(anndata)

cat("Run command\n")
input <- paste0(
  meta$resources_dir,
  "/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5"
)
output <- "out.h5"

system(paste0(
  "./", meta$name, " ",
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
expect_gte(out_ad$n_obs, 100)
expect_gte(out_ad$n_vars, 100)

# check if filtering has been performed
expect_lte(out_ad$n_obs, orig_ad$n_obs)
expect_lte(out_ad$n_vars, orig_ad$n_vars)

cat(">> All tests passed successfully!\n")
