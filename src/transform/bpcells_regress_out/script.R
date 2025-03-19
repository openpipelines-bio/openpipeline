cat("Loading libraries\n")
library(glue)
library(BPCells)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
mudata <- reticulate::import("mudata")

## VIASH START
par <- list(
  input = paste0(
    "resources_test/pbmc_1k_protein_v3",
    "pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
  ),
  output = "output.h5mu",
  modality = "rna"
)
## VIASH END

# Read the h5mu file and make var names unique
mdata <- mudata$read_h5mu(par$input)

# Regress out
if (!is.null(par$obs_keys) && length(par$obs_keys) > 0) {
  cat("Regress out variables ", par$obs_keys, " on modality ",
    par$modality, "\n",
    sep = ""
  )

  # Fetch modality AnnData and convert to an iterable matrix
  adata <- mdata$mod[[par$modality]]

  # Fetch the input layer
  mat <-
    if (is.null(par$input_layer)) {
      cat("Using .X as input layer\n")
      adata$X
    } else {
      cat("Using .layers ", par$input_layer, " as input layer\n", sep = "")
      adata$layers[[par$input_layer]]
    }

  imat <- as(as(mat, "CsparseMatrix"), "IterableMatrix")
  dimnames(imat) <- NULL

  # obs_keys is not NULL and not empty
  latent_data <- as.data.frame(adata$obs[, par$obs_keys])

  # Regress out using BPCells
  regressed_data <- regress_out(imat, latent_data, prediction_axis = "col")

  # Convert iterable matrix back to C sparse matrix
  rmat <- as(regressed_data, "dgCMatrix")

  # Assign regressed out data back to AnnData object
  if (is.null(par$output_layer)) {
    cat("Using .X as output layer\n")
    adata$X <- rmat
  } else {
    cat("Using .layers ", par$output_layer, " as output layer\n", sep = "")
    adata$layers[[par$output_layer]] <- rmat
  }
} else {
  cat("No obs_keys provided, skipping regression\n")
}

# Write to output h5mu file
mdata$write(par$output, compression = par$output_compression)
