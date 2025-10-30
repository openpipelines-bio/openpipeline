library(SingleR)
library(Matrix)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
mudata <- reticulate::import("mudata")

### VIASH START
par <- list(
  input = "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
  modality = "rna",
  input_layer = NULL,
  input_var_gene_names = "gene_symbol",
  input_reference_gene_overlap = 100,
  input_obs_clusters = NULL,
  reference = "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
  reference_layer = NULL,
  reference_var_input = NULL,
  reference_var_gene_names = NULL,
  # reference_var_gene_names = "ensemblid",
  reference_obs_target = "cell_ontology_class",
  output = "singler_output.h5mu",
  output_compression = "gzip",
  output_obs_predictions = "singler_labels",
  output_obs_probability = "singlr_proba",
  output_obsm_scores = "single_r_scores",
  output_obs_delta_next = "singler_delta_next",
  output_obs_pruned_predictions = "singler_pruned_labels",
  quantile = 0.8,
  fine_tune = TRUE,
  fine_tuning_threshold = 0.05,
  prune = TRUE,
  de_n_genes = NULL,
  de_method = "classic"
)
meta <- list(
  cpus = 4
)

### VIASH END

get_layer <- function(adata, layer, var_gene_names) {
  # find data
  data <-
    if (is.null(layer)) {
      adata$X
    } else {
      adata$layers[[layer]]
    }

  # check if data is available
  if (is.null(data)) {
    if (is.null(layer)) {
      stop("No layer specified and no .X slot available in the AnnData object.")
    } else {
      stop(
        "Requested layer '",
        layer,
        "' is not available in the AnnData object. Available layers: ",
        paste(names(adata$layers), collapse = ", ")
      )
    }
  }

  # Set matrix dimnames
  input_gene_names <- sanitize_ensembl_ids(adata, var_gene_names)
  dimnames(data) <- list(adata$obs_names, input_gene_names)

  # return output
  data
}

sanitize_ensembl_ids <- function(adata, gene_symbol = NULL) {
  if (is.null(gene_symbol)) {
    gene_names <- adata$var_names
  } else {
    gene_names <- adata$var[[gene_symbol]]
  }

  # Pattern matches Ensembl IDs: starts with ENS, followed by any characters,
  # then an eleven digit number, optionally followed by .version_number
  ensembl_pattern <- "^(ENS.*\\d{11})(?:\\.\\d+)?$"

  # Remove version numbers for ensembl ids only
  sanitized <- ifelse(
    grepl(ensembl_pattern, gene_names, perl = TRUE),
    gsub(ensembl_pattern, "\\1", gene_names, perl = TRUE),
    as.character(gene_names)
  )

  sanitized
}

subset_vars <- function(adata, subset_col) {
  # Check if column exists
  if (!subset_col %in% colnames(adata$var)) {
    stop(
      "Requested to use .var column '",
      subset_col,
      "' as a selection of genes, but the column is not available."
    )
  }

  # Get the column
  subset_values <- adata$var[[subset_col]]

  # Check for NA values
  if (any(is.na(subset_values))) {
    stop(
      "The .var column `",
      subset_col,
      "` contains NA values. Can not subset data."
    )
  }

  # Check if it's logical/boolean
  if (!is.logical(subset_values)) {
    stop(
      "Expected dtype of .var column '",
      subset_col,
      "' to be `logical`, but found ",
      class(subset_values)[1],
      ". Can not subset data."
    )
  }

  # Subset and return copy
  adata_subset <- adata[, subset_values]$copy()
  adata_subset
}

# Read input data
cat("Reading input file\n")
input_mdata <- mudata$read_h5mu(par$input)
input_adata <- input_mdata$mod[[par$modality]]
input_matrix <- Matrix::t(
  get_layer(input_adata, par$input_layer, par$input_var_gene_names)
)

# Read reference
cat("Reading reference file\n")
ref_adata <- mudata$read_h5ad(par$reference, mod = "rna")
if (!is.null(par$reference_var_input)) {
  ref_adata <- subset_vars(ref_adata, par$reference_var_input)
}
ref_matrix <- Matrix::t(
  get_layer(ref_adata, par$reference_layer, par$reference_var_gene_names)
)

# Check overlap genes
cat("Checking overlap of genes between input and reference file\n")
common_ens_ids <- intersect(rownames(ref_matrix), rownames(input_matrix))
if (length(common_ens_ids) < par$input_reference_gene_overlap) {
  stop(
    "The intersection of genes between the query and reference is too small.\n",
    paste0("Expected overlap: ", par$input_reference_gene_overlap),
    paste0("Detected overlap: ", length(common_ens_ids))
  )
}

# Calculate CPU cores
n_workers <- meta$cpus %||% max(1, parallel::detectCores() - 1)

# Assign target labels
if (!par$reference_obs_target %in% colnames(ref_adata$obs)) {
  stop(
    "Requested to use .obs column '",
    par$reference_obs_target,
    "' as target labels, but the column is not available."
  )
} else {
  target_labels <- ref_adata$obs[[par$reference_obs_target]]
}

cat("Performing SingleR cell type prediction\n")
predictions <- SingleR(
  test = input_matrix,
  ref = ref_matrix,
  labels = target_labels,
  clusters = par$input_obs_clusters,
  genes = "de",
  de.method = par$de_method,
  de.n = par$de_n_genes,
  quantile = par$quantile,
  fine.tune = par$fine_tune,
  tune.thresh = par$fine_tuning_threshold,
  prune = par$prune,
  num.threads = n_workers
)

cat("Writing output data\n")
# Writing output slots
input_adata$obs[[par$output_obs_predictions]] <-
  predictions$labels
input_adata$obs[[par$output_obs_probability]] <-
  apply(predictions$scores, 1, max)
input_adata$obs[[par$output_obs_delta_next]] <-
  predictions$delta.next
input_adata$obs[[par$output_obs_pruned_predictions]] <-
  predictions$pruned.labels
input_adata$obsm[[par$output_obsm_scores]] <-
  predictions$scores
# Writing output H5MU
input_mdata$write(par$output, compression = par$output_compression)
