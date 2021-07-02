# assume tidyverse is installed
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(testthat, warn.conflicts = FALSE)
library(Seurat)

cat("Checking whether output is correct\n")
out <- processx::run("./h5adtoseurat", c(
  "--input", "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",
  "--output", "output.rds"
))
  
expect_equal(out$status, 0)
expect_true(file.exists("output.rds"))
  
seurat <- readRDS(file = "output.rds")
  
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

