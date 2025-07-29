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

get_input_layer <- function(adata, input_layer = NULL) {
  if (is.null(input_layer)) {
    return(adata$X)
  } else {
    return(adata$layers[[input_layer]])
  }
}

sanitize_gene_names <- function(adata, gene_symbol = NULL) {
  if (is.null(gene_symbol)) {
    gene_names <- adata$var_names
  } else {
    gene_names <- adata$var[[gene_symbol]]
  }
  # Remove version numbers (dot followed by digits at end of string)
  sanitized <- gsub("\\.[0-9]+$", "", gene_names)
  return(sanitized)
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
  return(adata[, subset_values]$copy())
}

# Read input data
cat("Reading input file\n")
input_mdata <- mudata$read_h5mu(par$input)
input_adata <- input_mdata$mod[[par$modality]]
input_matrix <- Matrix::t(get_input_layer(input_adata, par$input_layer))
input_gene_names <- sanitize_gene_names(input_adata, par$input_var_gene_names)
dimnames(input_matrix) <- list(input_gene_names, input_adata$obs_names)

# Read reference
cat("Reading reference file\n")
ref_adata <- mudata$read_h5ad(par$reference, mod = "rna")
if (!is.null(par$reference_var_input)) {
  ref_adata <- subset_vars(ref_adata, par$reference_var_input)
}
ref_matrix <- Matrix::t(get_input_layer(ref_adata, par$reference_layer))
ref_gene_names <- sanitize_gene_names(ref_adata, par$reference_var_gene_names)
dimnames(ref_matrix) <- list(ref_gene_names, ref_adata$obs_names)

# Check overlap genes
cat("Checking overlap of genes between input and reference file\n")
common_ens_ids <- intersect(ref_gene_names, input_gene_names)
if (length(common_ens_ids) < par$input_reference_gene_overlap) {
  stop(
    "The intersection of genes between the query and reference is too small.\n",
    paste0("Expected overlap: ", par$input_reference_gene_overlap),
    paste0("Detected overlap: ", length(common_ens_ids))
  )
}

# Calculate CPU cores
if (!is.null(meta$cpus)) {
  n_workers <- meta$cpus
} else {
  n_workers <- max(1, parallel::detectCores() - 1)
}

cat("Performing SingleR cell type prediction\n")
predictions <- SingleR(
  test = input_matrix,
  ref = ref_matrix,
  labels = ref_adata$obs[[par$reference_obs_target]],
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
input_adata$obs[[par$output_obs_predictions]] <- predictions$labels
cat("pred success")
input_adata$obs[[par$output_obs_probability]] <- apply(predictions$scores, 1, max)
cat("proba success")
input_adata$obs[[par$output_obs_delta_next]] <- predictions$delta.next
cat("delta success")
input_adata$obs[[par$output_obs_pruned_predictions]] <- predictions$pruned.labels
cat("prune success")
input_adata$obsm[[par$output_obsm_scores]] <- predictions$scores
cat("scores success")
# Writing output H5MU
input_mdata$write(par$output, compression = par$output_compression)
