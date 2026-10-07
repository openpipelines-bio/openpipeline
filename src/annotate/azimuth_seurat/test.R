library(testthat, warn.conflicts = FALSE)
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))

## VIASH START
meta <- list(
  executable = "target/executable/annotate/azimuth_seurat/azimuth_seurat",
  resources_dir = "resources_test/pbmc_1k_protein_v3",
  name = "azimuth_seurat"
)
## VIASH END

input_file <- file.path(meta[["resources_dir"]], "pbmc_1k_protein_v3_mms.rds")
input <- readRDS(input_file)

run_component <- function(args) {
  processx::run(
    meta[["executable"]],
    args,
    error_on_status = FALSE,
    stderr_to_stdout = TRUE
  )
}

# Functions are namespaced here, as lintr's object_usage_linter checks
# function bodies and cannot resolve package functions in CI.
expect_input_preserved <- function(output) {
  testthat::expect_s4_class(output, "Seurat")
  testthat::expect_equal(
    SeuratObject::Cells(output),
    SeuratObject::Cells(input)
  )
  testthat::expect_equal(
    SeuratObject::DefaultAssay(output),
    SeuratObject::DefaultAssay(input)
  )
  testthat::expect_equal(
    SeuratObject::LayerData(output, assay = "RNA", layer = "X"),
    SeuratObject::LayerData(input, assay = "RNA", layer = "X")
  )
  testthat::expect_equal(
    SeuratObject::Layers(output[["RNA"]]),
    SeuratObject::Layers(input[["RNA"]])
  )
  testthat::expect_equal(output[["RNA"]][[]], input[["RNA"]][[]])
  testthat::expect_equal(output[[]][, colnames(input[[]])], input[[]])
  testthat::expect_true(all(
    SeuratObject::Reductions(input) %in% SeuratObject::Reductions(output)
  ))
  testthat::expect_equal(
    SeuratObject::Graphs(output),
    SeuratObject::Graphs(input)
  )
  testthat::expect_equal(names(output@misc), names(input@misc))
}

test_that("Annotation using gene symbols from the feature metadata", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--input_var_gene_names", "gene_symbol",
    "--reference", "pbmcref",
    "--output", output_file
  ))
  expect_equal(out$status, 0, info = out$stdout)
  expect_true(file.exists(output_file))

  output <- readRDS(output_file)
  expect_input_preserved(output)

  predicted_cols <- grep("^predicted\\.", colnames(output[[]]), value = TRUE)
  expected_cols <- c(
    "predicted.celltype.l1", "predicted.celltype.l1.score",
    "predicted.celltype.l2", "predicted.celltype.l3"
  )
  expect_true(all(expected_cols %in% predicted_cols))
  expect_false(any(is.na(output$predicted.celltype.l1)))
  # pbmcref annotations contain at least B, T and mono cell types
  expect_gt(length(unique(output$predicted.celltype.l1)), 3)

  expect_true("mapping.score" %in% colnames(output[[]]))
  expect_true(all(output$mapping.score >= 0 & output$mapping.score <= 1))

  expect_true(all(c("ref.umap", "integrated_dr") %in% Reductions(output)))
  umap <- Embeddings(output, "ref.umap")
  expect_equal(dim(umap), c(ncol(output), 2))
  expect_equal(rownames(umap), Cells(output))

  expect_true("prediction.score.celltype.l1" %in% Assays(output))
  scores <- LayerData(
    output,
    assay = "prediction.score.celltype.l1", layer = "data"
  )
  expect_equal(colnames(scores), Cells(output))
})

test_that("Annotation using Ensembl IDs with a subset of annotation levels", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--reference", "pbmcref",
    "--annotation_levels", "celltype.l1",
    "--do_adt",
    "--umap_name", "azimuth_umap",
    "--output_compression", "none",
    "--output", output_file
  ))
  expect_equal(out$status, 0, info = out$stdout)
  expect_true(file.exists(output_file))

  output <- readRDS(output_file)
  expect_input_preserved(output)

  expect_true("predicted.celltype.l1" %in% colnames(output[[]]))
  expect_false("predicted.celltype.l2" %in% colnames(output[[]]))
  expect_true("azimuth_umap" %in% Reductions(output))
  expect_false("ref.umap" %in% Reductions(output))
  expect_true("impADT" %in% Assays(output))
  expect_false("prediction.score.celltype.l2" %in% Assays(output))
})

test_that("Fails on a missing input layer", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--reference", "pbmcref",
    "--output", output_file
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "Layer 'counts' not found in assay 'RNA'")
  expect_false(file.exists(output_file))
})

test_that("Fails on an invalid reference", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--reference", "not_a_real_reference",
    "--output", output_file
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "Could not find a reference for")
  expect_false(file.exists(output_file))
})
