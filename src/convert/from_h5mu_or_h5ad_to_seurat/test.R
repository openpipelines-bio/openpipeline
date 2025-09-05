library(testthat, warn.conflicts = FALSE)
library(hdf5r)


## VIASH START
meta <- list(
  executable = "target/executable/convert/from_h5ad_to_seurat",
  resources_dir = "resources_test",
  name = "from_h5ad_to_seurat"
)
## VIASH END

cat("> Checking conversion of h5ad file\n")

in_h5ad <- paste0(
  meta[["resources_dir"]],
  "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna.h5ad"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", in_h5ad,
    "--output", out_rds
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)

cat("> Checking whether Seurat object is in the right format\n")
expect_is(obj, "Seurat")
expect_equal(names(slot(obj, "assays")), "RNA")

open_file <- H5File$new(in_h5ad, mode = "r+")

dim_rds <- dim(obj)
dim_ad <- open_file[["X"]]$attr_open("shape")$read()

expect_equal(dim_rds[1], dim_ad[2])
expect_equal(dim_rds[2], dim_ad[1])

cat("> Checking conversion of h5mu file\n")

in_h5mu <- paste0(
  meta[["resources_dir"]],
  "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)
out_rds <- "output.rds"

cat("> Running ", meta[["name"]], "\n", sep = "")
out <- processx::run(
  meta[["executable"]],
  c(
    "--input", in_h5mu,
    "--output", out_rds
  )
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))

cat("> Reading output file\n")
obj <- readRDS(file = out_rds)

cat("> Checking whether Seurat object is in the right format\n")
expect_is(obj, "Seurat")
expect_equal(names(slot(obj, "assays")), "RNA")

dim_rds <- dim(obj)
mu_in <- H5File$new(in_h5mu, mode = "r")
dim_ad <- mu_in[["/mod/rna/X"]]$attr_open("shape")$read()

expect_equal(dim_rds[1], dim_ad[2])
expect_equal(dim_rds[2], dim_ad[1])
