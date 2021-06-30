library(testthat)
requireNamespace("Seurat")

system(paste0(
  "./download_10x_dataset ",
  "--input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 ",
  "--output pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 ",
  "--min_library_size 1000 ",
  "--min_cells_per_gene 300"
))

seu <- Seurat::Read10X_h5("pbmc_1k_protein_v3_raw_feature_bc_matrix.h5")

gex <- seu$`Gene Expression`
expect_is(gex, "Matrix")
expect_gte(nrow(gex), 100)
expect_gte(ncol(gex), 100)
expect_gte(sum(gex), 100)

cat(">> All tests passed successfully!\n")
