library(testthat, warn.conflicts = FALSE)

## VIASH START
meta <- list(
  executable = "target/docker/convert/from_h5mu_to_seurat/from_h5mu_to_seurat",
  resources_dir = "resources_test",
  name = "from_h5mu_to_seurat"
)
## VIASH END

cat("> Checking whether output is correct\n")

in_h5mu <- paste0(
  meta[["resources_dir"]],
  "/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_mms.h5mu"
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
expect_equal(sort(names(slot(obj, "assays"))), sort(c("prot", "rna")))

obj_rna <- slot(obj, "assays")$rna
obj_prot <- slot(obj, "assays")$prot

# todo: check whether obj_rna and obj_prot have correct properties
