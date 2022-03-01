library(testthat)
requireNamespace("hdf5r", quietly = TRUE)

cat("Run command\n")
system(paste0(
  "./", meta$functionality_name, " ",
  "--input '", meta$resources_dir, "/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5' ",
  "--output out.h5 ",
  "--min_library_size 1000 ",
  "--min_cells_per_gene 300 ",
  "--verbose"
))

orig_seu <- Seurat::Read10X_h5("pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5")
orig_gex <- orig_seu$`Gene Expression`

seu <- Seurat::Read10X_h5("out.h5")
gex <- seu$`Gene Expression`

# check if matrix is not empty
expect_is(gex, "Matrix")
expect_gte(nrow(gex), 100)
expect_gte(ncol(gex), 100)
expect_gte(sum(gex), 100)

# check if filtering has been performed
expect_lte(nrow(gex), nrow(orig_gex))
expect_lte(ncol(gex), ncol(orig_gex))
expect_lte(sum(gex), sum(orig_gex))

cat(">> All tests passed successfully!\n")
