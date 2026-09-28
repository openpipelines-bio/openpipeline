cat("Loading libraries\n")
library(glue)
library(BPCells)
requireNamespace("reticulate", quietly = TRUE)
h5py <- reticulate::import("h5py")
anndata <- reticulate::import("anndata")

## VIASH START
par <- list(
  input = paste0(
    "resources_test/pbmc_1k_protein_v3/",
    "pbmc_1k_protein_v3_mms.h5mu"
  ),
  output = "output.h5mu",
  modality = "rna",
  obs_keys = "total_counts",
  var_input = NULL
)
## VIASH END

# Start from a copy of the input, the regressed data is written into it
invisible(file.copy(par$input, par$output, overwrite = TRUE))

# Regress out
if (!is.null(par$obs_keys) && length(par$obs_keys) > 0) {
  cat("Regress out variables ", par$obs_keys, " on modality ",
    par$modality, "\n",
    sep = ""
  )

  mod_path <- file.path("mod", par$modality)

  # Fetch the input layer
  input_path <-
    if (is.null(par$input_layer)) {
      cat("Using .X as input layer\n")
      file.path(mod_path, "X")
    } else {
      cat("Using .layers ", par$input_layer, " as input layer\n", sep = "")
      file.path(mod_path, "layers", par$input_layer)
    }

  # Read the input layer lazily from disk (genes x cells)
  imat <- open_matrix_anndata_hdf5(par$input, group = input_path)
  dimnames(imat) <- NULL
  cat("Input matrix: ", ncol(imat), " cells x ", nrow(imat), " genes\n",
    sep = ""
  )

  # Only read obs and var from the modality
  h5 <- h5py$File(par$input, "r")
  obs <- anndata$io$read_elem(h5[[file.path(mod_path, "obs")]])
  var <- anndata$io$read_elem(h5[[file.path(mod_path, "var")]])
  h5$close()

  # obs_keys is not NULL and not empty
  latent_data <- as.data.frame(obs[, par$obs_keys, drop = FALSE])
  # regress_out builds a formula from the column names
  colnames(latent_data) <- make.names(colnames(latent_data), unique = TRUE)

  mask_var <- NULL
  if (!is.null(par$var_input)) {
    mask_var <- as.logical(var[[par$var_input]])
    cat("Regressing out on ", sum(mask_var), " genes selected by .var column ",
      par$var_input, "\n",
      sep = ""
    )
    imat <- imat[mask_var, ]
  }

  # Regress out using BPCells
  regressed_data <- regress_out(imat, latent_data, prediction_axis = "row")

  # Non-selected genes are set to 0
  if (!is.null(mask_var)) {
    zeros <- as(
      Matrix::sparseMatrix(
        i = integer(0), j = integer(0), x = numeric(0),
        dims = c(sum(!mask_var), ncol(regressed_data))
      ),
      "IterableMatrix"
    )
    regressed_data <- rbind(regressed_data, zeros)
    regressed_data <- regressed_data[
      order(c(which(mask_var), which(!mask_var))),
    ]
  }

  output_path <-
    if (is.null(par$output_layer)) {
      cat("Using .X as output layer\n")
      file.path(mod_path, "X")
    } else {
      cat("Using .layers ", par$output_layer, " as output layer\n", sep = "")
      file.path(mod_path, "layers", par$output_layer)
    }

  # Remove the layer that will be overwritten
  h5 <- h5py$File(par$output, "r+")
  if (h5$`__contains__`(output_path)) h5$`__delitem__`(output_path)
  h5$close()

  # BPCells only supports gzip compression
  gzip_level <- if (is.null(par$output_compression)) 0L else 4L
  cat("Writing ", ncol(regressed_data), " cells x ", nrow(regressed_data),
    " genes with ",
    format(as.numeric(ncol(regressed_data)) * nrow(imat),
      big.mark = ",", scientific = FALSE
    ),
    " stored values, gzip level ", gzip_level, "\n",
    sep = ""
  )

  # Regressed values are computed while they are written to disk
  invisible(write_matrix_anndata_hdf5(
    regressed_data, par$output,
    group = output_path, gzip_level = gzip_level
  ))
} else {
  cat("No obs_keys provided, skipping regression\n")
}
