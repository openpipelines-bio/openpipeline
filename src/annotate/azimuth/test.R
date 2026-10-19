library(testthat, warn.conflicts = FALSE)
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))

## VIASH START
meta <- list(
  executable = "target/executable/annotate/azimuth/azimuth",
  resources_dir = "resources_test/pbmc_1k_protein_v3",
  name = "azimuth"
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

  # Without --annotation_levels, all levels of pbmcref are transferred and
  # stored under the default output names
  expected_cols <- c(
    paste0("azimuth_pred_celltype_", c("l1", "l2", "l3")),
    paste0("azimuth_probability_celltype_", c("l1", "l2", "l3")),
    "azimuth_mapping_score"
  )
  expect_true(all(expected_cols %in% colnames(output[[]])))
  expect_false(any(grepl("^predicted\\.", colnames(output[[]]))))
  expect_false("mapping.score" %in% colnames(output[[]]))
  expect_false(any(is.na(output$azimuth_pred_celltype_l1)))
  # pbmcref annotations contain at least B, T and mono cell types
  expect_gt(length(unique(output$azimuth_pred_celltype_l1)), 3)
  # Prediction scores can exceed 1 by a floating point rounding error
  expect_true(all(
    output$azimuth_probability_celltype_l1 >= 0 &
      output$azimuth_probability_celltype_l1 <= 1 + 1e-8
  ))
  expect_true(all(
    output$azimuth_mapping_score >= 0 & output$azimuth_mapping_score <= 1
  ))

  expect_true(all(
    c("X_azimuth_umap", "X_integrated_azimuth") %in% Reductions(output)
  ))
  expect_false(any(c("ref.umap", "integrated_dr") %in% Reductions(output)))
  umap <- Embeddings(output, "X_azimuth_umap")
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
    "--output_obs_predictions", "custom_pred",
    "--output_obs_probability", "custom_probability",
    "--output_obs_mapping_score", "custom_mapping_score",
    "--output_obsm_umap", "custom_umap",
    "--output_obsm_integrated", "custom_integrated",
    "--output_compression", "none",
    "--output", output_file
  ))
  expect_equal(out$status, 0, info = out$stdout)
  expect_true(file.exists(output_file))

  output <- readRDS(output_file)
  expect_input_preserved(output)

  expected_cols <- c(
    "custom_pred", "custom_probability", "custom_mapping_score"
  )
  expect_true(all(expected_cols %in% colnames(output[[]])))
  expect_false(any(grepl("^azimuth_", colnames(output[[]]))))
  expect_false(any(grepl("^predicted\\.", colnames(output[[]]))))
  expect_true(all(
    c("custom_umap", "custom_integrated") %in% Reductions(output)
  ))
  expect_false(any(
    c("X_azimuth_umap", "ref.umap", "integrated_dr") %in% Reductions(output)
  ))
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
  # --reference is restricted to a list of choices
  expect_match(out$stdout, "not_a_real_reference")
  expect_match(out$stdout, "choices")
  expect_false(file.exists(output_file))
})

test_that("Fails on --reference custom without --custom_reference", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--reference", "custom",
    "--output", output_file
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "--custom_reference must be provided")
  expect_false(file.exists(output_file))
})

test_that("Fails on --custom_reference without --reference custom", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--reference", "pbmcref",
    "--custom_reference", tempdir(),
    "--output", output_file
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "can only be used with --reference 'custom'")
  expect_false(file.exists(output_file))
})

test_that("Fails on a custom reference directory without reference files", {
  output_file <- tempfile(fileext = ".rds")
  empty_dir <- tempfile()
  dir.create(empty_dir)
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--reference", "custom",
    "--custom_reference", empty_dir,
    "--output", output_file
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "missing the required file")
  expect_false(file.exists(output_file))
})

test_that("Fails on output names without annotation levels", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--reference", "pbmcref",
    "--output_obs_predictions", "custom_pred",
    "--output", output_file
  ))
  expect_false(out$status == 0)
  expect_match(
    out$stdout,
    "--output_obs_predictions can only be used together with"
  )
  expect_false(file.exists(output_file))
})

test_that("Fails on output names not matching the annotation levels", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--reference", "pbmcref",
    "--annotation_levels", "celltype.l1;celltype.l2",
    "--output_obs_probability", "custom_probability",
    "--output", output_file
  ))
  expect_false(out$status == 0)
  expect_match(
    out$stdout,
    "--output_obs_probability should have the same length as"
  )
  expect_false(file.exists(output_file))
})

test_that("Fails on identical reduction names", {
  output_file <- tempfile(fileext = ".rds")
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--reference", "pbmcref",
    "--output_obsm_umap", "X_azimuth",
    "--output_obsm_integrated", "X_azimuth",
    "--output", output_file
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "should be different")
  expect_false(file.exists(output_file))
})
