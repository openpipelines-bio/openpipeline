library(decontX)
library(rhdf5)

### VIASH START
par <- list(
  input = "",
  modality = "rna",
  input_layer = NULL,
  input_obs_clusters = NULL,
  input_obs_batch = NULL,
  background = NULL,
  background_layer = NULL,
  background_obs_batch = NULL,
  max_iter = 500,
  delta = c(10, 10),
  estimate_delta = TRUE,
  convergence = 0.001,
  iter_log_lik = 10,
  var_genes = 5000,
  dbscan_eps = 1,
  seed = 12345,
  logfile = NULL,
  verbose = TRUE,
  output = "",
  output_layer = "decontx_counts",
  output_obs_contamination = "decontx_contamination",
  output_obs_clusters = "decontx_clusters",
  output_compression = NULL
)

meta <- list(
  cpus = 4,
  resources_dir = ".",
  temp_dir = tempdir()
)

### VIASH END

# Each modality of an .h5mu file is an AnnData group under /mod/<modality>.
# The input is read with anndataR after extracting the modality to a
# temporary .h5ad file. The output is written in place with rhdf5, following
# the AnnData on-disk format, see
# https://anndata.readthedocs.io/en/latest/fileformat-prose.html

# Uses rhdf5 rather than hdf5r: anndataR reads through rhdf5, and a file left
# open by a different HDF5 library stays locked for it.
h5mu_to_h5ad <- function(h5mu_path, modality_name) {
  tmp_path <- tempfile(fileext = ".h5ad")
  mod_path <- paste0("/mod/", modality_name)
  src <- H5Fopen(h5mu_path, flags = "H5F_ACC_RDONLY")
  on.exit(H5Fclose(src))
  if (!H5Lexists(src, mod_path)) {
    stop(
      "Modality '", modality_name, "' is not available in '", h5mu_path, "'."
    )
  }
  dest <- H5Fcreate(tmp_path)
  on.exit(H5Fclose(dest), add = TRUE)
  mod <- H5Gopen(src, mod_path)
  on.exit(H5Gclose(mod), add = TRUE, after = FALSE)

  for (child in h5ls(mod, recursive = FALSE)$name) {
    H5Ocopy(mod, child, dest, child)
  }
  # The modality group's attributes (encoding-type: anndata, ...) become the
  # root attributes of the .h5ad file
  h5_write_attrs(dest, "/", h5readAttributes(src, mod_path))

  tmp_path
}

# Read the counts matrix (as a genes x cells dgCMatrix) and the requested
# .obs columns of one modality
read_modality <- function(file, modality, layer, obs_columns = list()) {
  h5ad_path <- h5mu_to_h5ad(file, modality)
  on.exit(unlink(h5ad_path))
  adata <- anndataR::read_h5ad(h5ad_path, as = "InMemoryAnnData")

  mat <- if (is.null(layer)) adata$X else adata$layers[[layer]]
  if (is.null(mat)) {
    what <- if (is.null(layer)) ".X" else paste0("layer '", layer, "'")
    stop("Counts matrix ", what, " is not available in '", file, "'.")
  }
  mat <- Matrix::t(methods::as(mat, "CsparseMatrix"))
  dimnames(mat) <- list(adata$var_names, adata$obs_names)

  list(
    matrix = mat,
    obs = lapply(obs_columns, function(col) {
      if (is.null(col)) {
        return(NULL)
      }
      if (!col %in% colnames(adata$obs)) {
        stop(".obs column '", col, "' is not available in '", file, "'.")
      }
      # Categoricals are read as factors, decontX expects plain vectors
      as.vector(adata$obs[[col]])
    })
  )
}

# Writers for the AnnData elements added to the output
h5_write_attrs <- function(fid, path, attrs, scalar = TRUE) {
  oid <- H5Oopen(fid, path)
  on.exit(H5Oclose(oid))
  for (name in names(attrs)) {
    if (H5Aexists(oid, name)) H5Adelete(oid, name)
    value <- attrs[[name]]
    if (is.character(value)) {
      h5writeAttribute(
        value, oid, name,
        encoding = "UTF-8", variableLengthString = TRUE, asScalar = scalar
      )
    } else {
      h5writeAttribute(value, oid, name, asScalar = scalar)
    }
  }
}

h5_delete_if_exists <- function(fid, path) {
  if (H5Lexists(fid, path)) h5delete(fid, path)
}

compression_type <- if (is.null(par$output_compression)) {
  "none"
} else {
  par$output_compression
}
compression <- switch(compression_type,
  gzip = list(filter = "GZIP", level = 4),
  lzf = list(filter = "LZF", level = 0),
  none = list(filter = "NONE", level = 0)
)

h5_write_array <- function(fid, path, value, encoding_type = "array") {
  n <- length(value)
  # Filters don't apply to variable-length strings (LZF even segfaults on them)
  filter <- if (is.character(value)) "NONE" else compression$filter
  h5createDataset(
    fid, path,
    dims = n,
    storage.mode = storage.mode(value),
    size = NULL,
    encoding = if (is.character(value)) "UTF-8" else NULL,
    chunk = max(1L, min(n, 262144L)),
    filter = filter,
    level = if (filter == "NONE") 0 else compression$level
  )
  h5write(value, fid, path)
  h5_write_attrs(
    fid, path,
    list("encoding-type" = encoding_type, "encoding-version" = "0.2.0")
  )
}

# mat is genes x cells dgCMatrix, stored as a cells x genes CSR matrix
h5_write_csr_matrix <- function(fid, path, mat) {
  mat <- methods::as(mat, "CsparseMatrix")
  h5_delete_if_exists(fid, path)
  h5createGroup(fid, path)
  h5_write_array(fid, paste0(path, "/data"), as.numeric(mat@x))
  h5_write_array(fid, paste0(path, "/indices"), mat@i)
  h5_write_array(fid, paste0(path, "/indptr"), mat@p)
  h5_write_attrs(
    fid, path,
    list("encoding-type" = "csr_matrix", "encoding-version" = "0.1.0")
  )
  h5_write_attrs(fid, path, list(shape = rev(dim(mat))), scalar = FALSE)
}

h5_write_categorical <- function(fid, path, values) {
  values <- factor(values)
  h5_delete_if_exists(fid, path)
  h5createGroup(fid, path)
  codes <- as.integer(values) - 1L
  codes[is.na(codes)] <- -1L
  h5_write_array(fid, paste0(path, "/codes"), codes)
  h5_write_array(
    fid, paste0(path, "/categories"), levels(values),
    encoding_type = "string-array"
  )
  h5_write_attrs(
    fid, path,
    list(
      "encoding-type" = "categorical",
      "encoding-version" = "0.2.0",
      ordered = FALSE
    )
  )
}

# Read input data
cat("Reading input file\n")
input <- read_modality(
  par$input, par$modality, par$input_layer,
  list(clusters = par$input_obs_clusters, batch = par$input_obs_batch)
)

# Read background data, if provided
background <- NULL
if (!is.null(par$background)) {
  cat("Reading background file\n")
  background_layer <- if (is.null(par$background_layer)) {
    par$input_layer
  } else {
    par$background_layer
  }
  background <- read_modality(
    par$background, par$modality, background_layer,
    list(batch = par$background_obs_batch)
  )
}

cat("Estimating and removing contamination with DecontX\n")
result <- decontX(
  x = input$matrix,
  z = input$obs$clusters,
  batch = input$obs$batch,
  background = background$matrix,
  bgBatch = background$obs$batch,
  maxIter = par$max_iter,
  delta = par$delta,
  estimateDelta = par$estimate_delta,
  convergence = par$convergence,
  iterLogLik = par$iter_log_lik,
  varGenes = par$var_genes,
  dbscanEps = par$dbscan_eps,
  seed = par$seed,
  logfile = par$logfile,
  verbose = par$verbose
)

cat("Writing output data\n")
# Copy the input and add the results to it in place, leaving every other
# modality and field untouched.
if (!file.copy(par$input, par$output, overwrite = TRUE)) {
  stop("Could not copy '", par$input, "' to '", par$output, "'.")
}
fid <- H5Fopen(par$output, flags = "H5F_ACC_RDWR")
mod_path <- paste0("/mod/", par$modality)
obs_path <- paste0(mod_path, "/obs")

h5_write_csr_matrix(
  fid, paste0(mod_path, "/layers/", par$output_layer), result$decontXcounts
)

contamination_path <- paste0(obs_path, "/", par$output_obs_contamination)
h5_delete_if_exists(fid, contamination_path)
h5_write_array(fid, contamination_path, as.numeric(result$contamination))
h5_write_categorical(
  fid, paste0(obs_path, "/", par$output_obs_clusters), result$z
)

# New .obs columns are only read when listed in the 'column-order' attribute
column_order <- as.character(h5readAttributes(fid, obs_path)[["column-order"]])
column_order <- union(
  column_order, c(par$output_obs_contamination, par$output_obs_clusters)
)
h5_write_attrs(
  fid, obs_path, list("column-order" = column_order),
  scalar = FALSE
)

H5Fclose(fid)
