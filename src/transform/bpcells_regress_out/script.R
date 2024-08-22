cat("Loading libraries\n")
library(glue)
library(BPCells)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("reticulate")
library(assertthat)
mudata <- reticulate::import("mudata")

## VIASH START
par <- list(
    input = "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    output = "output.h5mu",
    modality = "rna"
)
## VIASH END

# Read the h5mu file and make var names unique
mdata <- mudata$read_h5mu(par$input)
mdata$var_names_make_unique()

# Fetch modality AnnData and convert to an iterable matrix
adata <- mdata$mod[[par$modality]]
mat <- as.matrix(adata)
imat <- as(as(mat, "CsparseMatrix"), "IterableMatrix")

# Regress out
if (!is.null(par$obs_keys) && length(par$obs_keys) > 0) {
  glue("Regress out variables {par$obs_keys} on modality {par$modality}")
  # obs_keys is not NULL and not empty
  latent_data <- as.data.frame(adata$obs[, par$obs_keys])
  # Regress out using BPCells
  regressed_data <- regress_out(imat, latent_data, prediction_axis = "col")
  # Convert iterable matrix back to R sparse matrix
  rmat <- as(regressed_data, "dgCMatrix")
}

# Assign regressed out data back to AnnData object
adata$X <- rmat
mdata$mod[[par$modality]] <- adata

# Write to output h5mu file
mdata$write(par$output, compression=par$output_compression)
