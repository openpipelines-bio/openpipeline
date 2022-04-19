# assume tidyverse is installed
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, warn.conflicts = FALSE)
library(Seurat)

cat("Checking whether output is correct\n")

in_h5mu <- paste0(meta_resources_dir, "/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tx_processing.h5mu")
in_h5ad <- "temp.h5ad"
out_rds <- "output.rds"

out <- processx::run(
  paste0("./", meta_functionality_name),
  c("--input", in_h5ad, "--output", out_rds)
)

expect_equal(out$status, 0)
expect_true(file.exists(out_rds))
  
seurat <- readRDS(file = out_rds)
  
  #class(seurat)
 # expect_is(seurat, "ExpressionSet")
  
 # pdata <- as(Biobase::phenoData(eset), "data.frame")
 # fdata <- as(Biobase::featureData(eset), "data.frame")
 # counts <- Biobase::assayData(eset)$counts
 # exprs <- Biobase::assayData(eset)$exprs
  
 # expect_equal(rownames(counts), rownames(exprs))
 # expect_equal(rownames(counts), rownames(fdata))
  
 # expect_equal(colnames(counts), colnames(exprs))
 # expect_equal(colnames(counts), rownames(pdata))

