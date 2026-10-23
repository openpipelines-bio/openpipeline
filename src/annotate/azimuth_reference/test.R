library(testthat, warn.conflicts = FALSE)
suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))

## VIASH START
meta <- list(
  executable = "target/executable/annotate/azimuth_reference/azimuth_reference",
  resources_dir = "resources_test/annotation_test_data",
  name = "azimuth_reference"
)
## VIASH END

input_file <- file.path(meta[["resources_dir"]], "TS_Blood_filtered.rds")
input <- readRDS(input_file)

run_component <- function(args) {
  processx::run(
    meta[["executable"]],
    args,
    error_on_status = FALSE,
    stderr_to_stdout = TRUE
  )
}

expect_valid_reference <- function(output_dir, annotation_levels) {
  testthat::expect_true(dir.exists(output_dir))
  testthat::expect_setequal(
    list.files(output_dir),
    c("ref.Rds", "idx.annoy")
  )
  # Load the reference the same way annotate/azimuth does
  ref <- Azimuth::LoadReference(output_dir)$map
  testthat::expect_s4_class(ref, "Seurat")
  testthat::expect_equal(SeuratObject::Cells(ref), SeuratObject::Cells(input))
  testthat::expect_true(
    all(c("refDR", "refUMAP") %in% SeuratObject::Reductions(ref))
  )
  testthat::expect_true("refAssay" %in% SeuratObject::Assays(ref))
  testthat::expect_true(all(annotation_levels %in% colnames(ref[[]])))
  for (level in annotation_levels) {
    testthat::expect_equal(
      as.character(ref[[level, drop = TRUE]]),
      as.character(input[[level, drop = TRUE]])
    )
  }
  ref
}

test_that("Build a reference with multiple annotation levels", {
  output_dir <- tempfile()
  out <- run_component(c(
    "--input", input_file,
    "--annotation_levels", "cell_ontology_class;cell_type",
    "--output", output_dir
  ))
  expect_equal(out$status, 0, info = out$stdout)

  ref <- expect_valid_reference(
    output_dir, c("cell_ontology_class", "cell_type")
  )
  expect_equal(Azimuth::ReferenceVersion(ref), "0.0.0")
  # The feature names of the input (gene symbols) are used, in which
  # Seurat replaces underscores with dashes
  expect_true(all(rownames(ref) %in% gsub("_", "-", rownames(input))))
  expect_equal(nrow(ref), 3000)
  neighbors <- ref[["refdr.annoy.neighbors"]]
  expect_equal(neighbors@alg.info$ndim, 50)
  expect_equal(dim(Indices(neighbors)), c(ncol(ref), 31))
  expect_equal(ncol(Embeddings(ref, "refDR")), 50)
})

test_that("Build a reference with custom reference building arguments", {
  output_dir <- tempfile()
  out <- run_component(c(
    "--input", input_file,
    "--annotation_levels", "cell_ontology_class",
    "--reference_version", "1.2.3",
    "--n_variable_features", "2000",
    "--n_dims", "60",
    "--k_param", "20",
    "--seed", "1",
    "--output", output_dir
  ))
  expect_equal(out$status, 0, info = out$stdout)

  ref <- expect_valid_reference(output_dir, "cell_ontology_class")
  expect_false("cell_type" %in% colnames(ref[[]]))
  expect_equal(Azimuth::ReferenceVersion(ref), "1.2.3")
  expect_equal(nrow(ref), 2000)
  neighbors <- ref[["refdr.annoy.neighbors"]]
  expect_equal(neighbors@alg.info$ndim, 60)
  expect_equal(dim(Indices(neighbors)), c(ncol(ref), 20))
  expect_equal(ncol(Embeddings(ref, "refDR")), 60)
})

test_that("Build a reference using gene names from the feature metadata", {
  output_dir <- tempfile()
  out <- run_component(c(
    "--input", input_file,
    "--input_var_gene_names", "ensemblid",
    "--annotation_levels", "cell_ontology_class",
    "--output", output_dir
  ))
  expect_equal(out$status, 0, info = out$stdout)

  ref <- expect_valid_reference(output_dir, "cell_ontology_class")
  # Features without an Ensembl ID keep their original name
  ensembl_ids <- input[["RNA"]][[]]$ensemblid
  expected_names <- ifelse(is.na(ensembl_ids), rownames(input), ensembl_ids)
  expect_true(all(rownames(ref) %in% gsub("_", "-", expected_names)))
  expect_gt(mean(grepl("^ENSG", rownames(ref))), 0.9)
})

test_that("Fails on a missing input layer", {
  output_dir <- tempfile()
  out <- run_component(c(
    "--input", input_file,
    "--input_layer", "X",
    "--annotation_levels", "cell_ontology_class",
    "--output", output_dir
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "Layer 'X' not found in assay 'RNA'")
  expect_false(dir.exists(output_dir))
})

test_that("Fails on a missing annotation level", {
  output_dir <- tempfile()
  out <- run_component(c(
    "--input", input_file,
    "--annotation_levels", "cell_ontology_class;not_a_column",
    "--output", output_dir
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "Annotation level\\(s\\) not found .*not_a_column")
  expect_false(dir.exists(output_dir))
})

test_that("Fails on an annotation level with missing labels", {
  output_dir <- tempfile()
  na_input <- input
  na_input$cell_ontology_class <- as.character(na_input$cell_ontology_class)
  na_input$cell_ontology_class[1:5] <- NA
  na_input_file <- tempfile(fileext = ".rds")
  saveRDS(na_input, na_input_file)

  out <- run_component(c(
    "--input", na_input_file,
    "--annotation_levels", "cell_ontology_class",
    "--output", output_dir
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "contains 5 cell\\(s\\) without a label")
  expect_false(dir.exists(output_dir))
})

test_that("Fails on more dimensions than cells", {
  output_dir <- tempfile()
  out <- run_component(c(
    "--input", input_file,
    "--annotation_levels", "cell_ontology_class",
    "--n_dims", "1000",
    "--output", output_dir
  ))
  expect_false(out$status == 0)
  expect_match(out$stdout, "should not be larger than the number of cells")
  expect_false(dir.exists(output_dir))
})
