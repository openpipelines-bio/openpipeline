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

cat("Using .layers ", par$output_layer, " as output layer\n", sep = "")
output_layer <- file.path(mod_path, "layers", par$output_layer)

# Read the input layer lazily from disk (genes x cells)
imat <- open_matrix_anndata_hdf5(par$input, group = input_layer)
dimnames(imat) <- NULL
cat("Input matrix: ", ncol(imat), " cells x ", nrow(imat), " genes\n",
  sep = ""
)

h5 <- h5py$File(par$input, "r")
# Only read obs and var from the modality
obs <- anndata$io$read_elem(h5[[file.path(mod_path, "obs")]])
var <- anndata$io$read_elem(h5[[file.path(mod_path, "var")]])
# Make sure output layer does not exist in the data yet
output_exists <- h5$`__contains__`(output_layer)
h5$close()
if (output_exists) {
  stop("Output layer ", par$output_layer, " already exists in modality ",
    par$modality, ", please choose a new layer name.",
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

# Regress out using BPCells
cat("Setting up regression with covariates: ",
  paste(colnames(obs), collapse = ", "), "\n",
  sep = ""
)
regressed_data <- regress_out(imat, obs, prediction_axis = "row")

# Non-selected genes are set to 0
if (!is.null(mask_var)) {
  cat("Setting ", sum(!mask_var), " non-selected genes to 0\n", sep = "")
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

# BPCells only supports gzip compression
gzip_level <-
  if (is.null(par$output_compression)) 0L else par$output_compression
cat("Writing regressed data to ", output_layer, " with gzip level ",
  gzip_level, "\n",
  sep = ""
)

# Regressed values are computed while they are written to disk
write_matrix_anndata_hdf5(
  regressed_data, par$output,
  group = output_layer, gzip_level = gzip_level
)
