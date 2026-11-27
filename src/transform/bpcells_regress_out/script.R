cat("Loading libraries\n")
library(glue)
library(BPCells)
requireNamespace("anndataR", quietly = TRUE)

# Reading a single modality from a MuData file via `group` and write_elem()
# require anndataR >= 1.4.0
stopifnot(packageVersion("anndataR") >= "1.4.0")

## VIASH START
par <- list(
  input = paste0(
    "resources_test/pbmc_1k_protein_v3/",
    "pbmc_1k_protein_v3_mms.h5mu"
  ),
  output = "output.h5mu",
  modality = "rna",
  obs_keys = "total_counts",
  var_input = NULL,
  output_layer = "regressed",
  obsm_pca_output = NULL,
  varm_pca_output = "regressed_pca_loadings",
  uns_pca_output = "regressed_pca_variance",
  num_components = 50L,
  scale_max_value = NULL,
  overwrite = FALSE
)
meta <- list(temp_dir = tempdir())
## VIASH END

# An empty --output_layer (e.g. `--output_layer ""` on the command line)
# disables storing the regressed data, the same as passing null from Nextflow
if (identical(par$output_layer, "")) {
  par$output_layer <- NULL
}
run_pca <- !is.null(par$obsm_pca_output)
if (is.null(par$output_layer) && !run_pca) {
  stop("At least one of --output_layer or --obsm_pca_output must be provided.",
    call. = FALSE
  )
}

# Start from a copy of the input, the regressed data is written directly to disk
file.copy(par$input, par$output, overwrite = TRUE)

# Regress out
cat("Regress out variables ", par$obs_keys, " on modality ",
  par$modality, "\n",
  sep = ""
)

mod_path <- file.path("mod", par$modality)

# Fetch the input layer
input_layer <-
  if (is.null(par$input_layer)) {
    cat("Using .X as input layer\n")
    file.path(mod_path, "X")
  } else {
    cat("Using .layers ", par$input_layer, " as input layer\n", sep = "")
    file.path(mod_path, "layers", par$input_layer)
  }

output_layer <- NULL
if (!is.null(par$output_layer)) {
  cat("Using .layers ", par$output_layer, " as output layer\n", sep = "")
  output_layer <- file.path(mod_path, "layers", par$output_layer)
}

# Read the input layer lazily from disk (genes x cells)
imat <- open_matrix_anndata_hdf5(par$input, group = input_layer)
dimnames(imat) <- NULL
cat("Input matrix: ", ncol(imat), " cells x ", nrow(imat), " genes\n",
  sep = ""
)

# Elements of the modality are only read from disk when they are accessed,
# so only obs and var are loaded
adata <- anndataR::read_h5ad(par$input, as = "HDF5AnnData", group = mod_path)
obs <- adata$obs
var <- adata$var
# Make sure output layer does not exist in the data yet. Listing the keys of
# a slot does not read its contents, and gives NULL when the slot is missing.
output_exists <- !is.null(par$output_layer) &&
  par$output_layer %in% adata$layers_keys()
# PCA slots can only be replaced when --overwrite is set
existing_pca_slots <- character(0)
if (run_pca) {
  pca_slots <- c(
    obsm = par$obsm_pca_output,
    varm = par$varm_pca_output,
    uns = par$uns_pca_output
  )
  slot_keys <- list(
    obsm = adata$obsm_keys(),
    varm = adata$varm_keys(),
    uns = adata$uns_keys()
  )
  pca_slot_exists <- vapply(
    names(pca_slots),
    function(slot) pca_slots[[slot]] %in% slot_keys[[slot]],
    logical(1)
  )
  existing_pca_slots <- file.path(
    mod_path, names(pca_slots), pca_slots
  )[pca_slot_exists]
}
if (output_exists) {
  stop("Output layer ", par$output_layer, " already exists in modality ",
    par$modality, ", please choose a new layer name.",
    call. = FALSE
  )
}
if (length(existing_pca_slots) > 0 && !isTRUE(par$overwrite)) {
  stop("PCA slot(s) ", paste(existing_pca_slots, collapse = ", "),
    " already exist, use --overwrite to replace them.",
    call. = FALSE
  )
}

# select and sanitize obs names to for regression formula
obs <- as.data.frame(obs[, par$obs_keys, drop = FALSE])
colnames(obs) <- make.names(colnames(obs), unique = TRUE)

# subset to HVG if requested
mask_var <- NULL
if (!is.null(par$var_input)) {
  mask_var <- as.logical(var[[par$var_input]])
  cat("Regressing out on ", sum(mask_var), " genes selected by .var column ",
    par$var_input, "\n",
    sep = ""
  )
  imat <- imat[mask_var, ]
} else {
  cat("No .var column provided, regressing out on all ", nrow(imat),
    " genes\n",
    sep = ""
  )
}

# The PCA reads the input many times, store it in the fast BPCells format
# instead of decompressing the h5 layer on every pass
if (run_pca) {
  cat("Storing input matrix in a temporary BPCells directory\n")
  imat <- write_matrix_dir(
    imat, tempfile("bpcells_input_", tmpdir = meta$temp_dir)
  )
}

# Regress out using BPCells
cat("Setting up regression with covariates: ",
  paste(colnames(obs), collapse = ", "), "\n",
  sep = ""
)
regressed_data <- regress_out(imat, obs, prediction_axis = "row")

if (!is.null(output_layer)) {
  output_data <- regressed_data
  # Non-selected genes are set to 0
  if (!is.null(mask_var)) {
    cat("Setting ", sum(!mask_var), " non-selected genes to 0\n", sep = "")
    zeros <- as(
      Matrix::sparseMatrix(
        i = integer(0), j = integer(0), x = numeric(0),
        dims = c(sum(!mask_var), ncol(output_data))
      ),
      "IterableMatrix"
    )
    output_data <- rbind(output_data, zeros)
    output_data <- output_data[
      order(c(which(mask_var), which(!mask_var))),
    ]
  }

  # BPCells only supports gzip compression
  gzip_level <-
    if (is.null(par$output_layer_compression)) {
      0L
    } else {
      par$output_layer_compression
    }
  cat("Writing regressed data to ", output_layer, " with gzip level ",
    gzip_level, "\n",
    sep = ""
  )

  # Regressed values are computed while they are written to disk
  write_matrix_anndata_hdf5(
    output_data, par$output,
    group = output_layer, gzip_level = gzip_level
  )
}

if (run_pca) {
  num_components <- par$num_components

  # Scale each gene to unit variance, genes without variance are only centered
  cat("Scaling regressed data\n")
  regressed_stats <- matrix_stats(regressed_data, row_stats = "variance")
  gene_mean <- regressed_stats$row_stats["mean", ]
  gene_variance <- regressed_stats$row_stats["variance", ]
  gene_sd <- sqrt(gene_variance)
  gene_sd[gene_sd == 0] <- 1
  scaled_data <- (regressed_data - gene_mean) / gene_sd
  # Every scaled gene with variance has unit variance
  total_variance <- sum(gene_variance > 0)

  if (!is.null(par$scale_max_value)) {
    max_value <- par$scale_max_value
    cat("Clipping scaled data to [-", max_value, ", ", max_value, "]\n",
      sep = ""
    )
    # BPCells only clips from above, the lower bound is applied to the
    # negated values
    scaled_data <- min_scalar(scaled_data, max_value)
    scaled_data <- min_scalar(scaled_data * -1, max_value) * -1

    # Clipping shifts the gene means and variances, so center again before
    # the PCA (as scanpy does)
    scaled_stats <- matrix_stats(scaled_data, row_stats = "variance")
    scaled_data <- scaled_data - scaled_stats$row_stats["mean", ]
    total_variance <- sum(scaled_stats$row_stats["variance", ])
  }

  # Matrix products with the lazy matrix are computed while streaming from disk
  cat("Computing ", num_components, " principal components\n", sep = "")
  pca <- svds(scaled_data, k = num_components)

  embedding <- sweep(pca$v, 2, pca$d, "*")
  loadings <- pca$u
  # Non-selected genes get zero loadings
  if (!is.null(mask_var)) {
    loadings <- matrix(0, nrow = length(mask_var), ncol = num_components)
    loadings[mask_var, ] <- pca$u
  }
  variance <- pca$d^2 / (ncol(scaled_data) - 1)

  cat("Writing PCA to .obsm ", par$obsm_pca_output,
    ", .varm ", par$varm_pca_output, " and .uns ", par$uns_pca_output, "\n",
    sep = ""
  )
  write_mod_elem <- function(slot, key, value) {
    # Existing slots were checked against --overwrite before the computation.
    # write_elem() creates missing parent groups.
    anndataR::write_elem(
      par$output, file.path(mod_path, slot, key), value,
      overwrite = TRUE
    )
  }
  write_mod_elem("obsm", par$obsm_pca_output, embedding)
  write_mod_elem("varm", par$varm_pca_output, loadings)
  write_mod_elem("uns", par$uns_pca_output, list(
    variance = variance,
    variance_ratio = variance / total_variance
  ))
}
