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

# The .h5mu files are read and written directly with rhdf5, following the
# AnnData on-disk format, see
# https://anndata.readthedocs.io/en/latest/fileformat-prose.html
# Each modality is an AnnData group under /mod/<modality>.

h5_encoding_type <- function(fid, path) {
  encoding <- h5readAttributes(fid, path)[["encoding-type"]]
  if (is.null(encoding)) "" else encoding
}

h5_check_exists <- function(fid, path, what) {
  if (!H5Lexists(fid, path)) {
    stop(what, " '", path, "' is not available in the input file.")
  }
}

# Read a counts matrix (.X or a layer) as a genes x cells dgCMatrix
read_matrix <- function(fid, path) {
  encoding <- h5_encoding_type(fid, path)
  if (encoding %in% c("csr_matrix", "csc_matrix")) {
    shape <- h5readAttributes(fid, path)[["shape"]]
    mat <- Matrix::sparseMatrix(
      i = as.integer(h5read(fid, paste0(path, "/indices"))),
      p = as.integer(h5read(fid, paste0(path, "/indptr"))),
      x = as.numeric(h5read(fid, paste0(path, "/data"))),
      # A cells x genes CSR matrix has the same layout as a genes x cells CSC
      dims = if (encoding == "csr_matrix") rev(shape) else shape,
      index1 = FALSE
    )
    if (encoding == "csr_matrix") mat else Matrix::t(mat)
  } else if (encoding == "array") {
    # rhdf5 reads the row-major cells x genes array as genes x cells
    methods::as(h5read(fid, path), "CsparseMatrix")
  } else {
    stop("Unsupported encoding '", encoding, "' for matrix '", path, "'.")
  }
}

# Read a single .obs column as a plain vector
read_obs_column <- function(fid, obs_path, column) {
  path <- paste0(obs_path, "/", column)
  h5_check_exists(fid, path, ".obs column")
  encoding <- h5_encoding_type(fid, path)
  if (encoding == "categorical") {
    codes <- as.integer(h5read(fid, paste0(path, "/codes")))
    categories <- as.vector(h5read(fid, paste0(path, "/categories")))
    codes[codes < 0] <- NA
    categories[codes + 1L]
  } else if (encoding %in% c("array", "string-array")) {
    as.vector(h5read(fid, path))
  } else {
    stop("Unsupported encoding '", encoding, "' for .obs column '", path, "'.")
  }
}

# Read the counts matrix and the requested .obs columns of one modality
read_modality <- function(file, modality, layer, obs_columns = list()) {
  fid <- H5Fopen(file, flags = "H5F_ACC_RDONLY")
  on.exit(H5Fclose(fid))

  mod_path <- paste0("/mod/", modality)
  h5_check_exists(fid, mod_path, "Modality")
  matrix_path <- if (is.null(layer)) {
    paste0(mod_path, "/X")
  } else {
    paste0(mod_path, "/layers/", layer)
  }
  h5_check_exists(fid, matrix_path, "Counts matrix")

  read_index <- function(df_path) {
    index_col <- h5readAttributes(fid, df_path)[["_index"]]
    as.character(h5read(fid, paste0(df_path, "/", index_col)))
  }
  obs_path <- paste0(mod_path, "/obs")

  mat <- read_matrix(fid, matrix_path)
  dimnames(mat) <- list(
    read_index(paste0(mod_path, "/var")),
    read_index(obs_path)
  )

  list(
    matrix = mat,
    obs = lapply(obs_columns, function(col) {
      if (is.null(col)) NULL else read_obs_column(fid, obs_path, col)
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
