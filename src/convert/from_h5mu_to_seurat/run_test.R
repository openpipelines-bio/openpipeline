library(testthat, warn.conflicts = FALSE)

cat("Checking whether output is correct\n")

# TODO: this should also work for ums and mms, but it doesn't.
in_h5mu <- paste0(meta[["resources_dir"]], "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_uss.h5mu")
out_rds <- "output.rds"

cat("> Running ", meta[["functionality_name"]], "\n", sep = "")
out <- processx::run(
  paste0("./", meta[["functionality_name"]]),
  c("--input", in_h5mu, "--output", out_rds)
)

cat("> Checking whether output file exists\n")
expect_equal(out$status, 0)
expect_true(file.exists(out_rds))
  
cat("> Reading output file\n")
obj <- readRDS(file = out_rds)

cat("> Checking whether Seurat object is in the right format\n")
expect_is(obj, "Seurat")
expect_equal(sort(names(slot(obj, "assays"))), sort(c("prot", "rna")))

obj_rna <- slot(obj, "assays")$rna
obj_prot <- slot(obj, "assays")$prot

# todo: check whether obj_rna and obj_prot have correct properties