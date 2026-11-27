library(decontX)

# Reading a single modality from a MuData file via `group` and the
# read_elem()/write_elem() functions require anndataR >= 1.4.0
stopifnot(packageVersion("anndataR") >= "1.4.0")

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
  verbose = FALSE,
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
# The modality is read directly from that group, and the results are written
# as single elements into a copy of the input, leaving every other modality
# and field untouched.

# Read the counts matrix (as a genes x cells dgCMatrix) and the requested
# .obs columns of one modality
read_modality <- function(file, modality, layer, obs_columns = list()) {
  group <- paste0("mod/", modality)
  # Elements are only read from disk when they are accessed
  adata <- anndataR::read_h5ad(file, as = "HDF5AnnData", group = group)

  # adata$layers would read all layers, only read the requested one
  mat <- if (is.null(layer)) {
    adata$X
  } else if (layer %in% adata$layers_keys()) {
    anndataR::read_elem(file, file.path(group, "layers", layer))
  }
  if (is.null(mat)) {
    what <- if (is.null(layer)) ".X" else paste0("layer '", layer, "'")
    stop("Counts matrix ", what, " is not available in '", file, "'.")
  }
  mat <- Matrix::t(methods::as(mat, "CsparseMatrix"))
  dimnames(mat) <- list(adata$var_names, adata$obs_names)

  obs <- adata$obs
  list(
    matrix = mat,
    obs = lapply(obs_columns, function(col) {
      if (is.null(col)) {
        return(NULL)
      }
      if (!col %in% colnames(obs)) {
        stop(".obs column '", col, "' is not available in '", file, "'.")
      }
      # Categoricals are read as factors, decontX expects plain vectors
      as.vector(obs[[col]])
    })
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
compression <- if (is.null(par$output_compression)) {
  "none"
} else {
  par$output_compression
}

# Copy the input and add the results to it in place, leaving every other
# modality and field untouched.
if (!file.copy(par$input, par$output, overwrite = TRUE)) {
  stop("Could not copy '", par$input, "' to '", par$output, "'.")
}
mod_path <- paste0("mod/", par$modality)

# write_elem() creates the layers group if it is missing
anndataR::write_elem(
  par$output,
  file.path(mod_path, "layers", par$output_layer),
  methods::as(Matrix::t(result$decontXcounts), "RsparseMatrix"),
  overwrite = TRUE,
  compression = compression
)

# Write the new .obs columns one by one: rewriting the whole .obs data frame
# could alter the dtypes of the existing columns
obs_path <- file.path(mod_path, "obs")
new_obs <- list(
  as.numeric(result$contamination),
  factor(result$z)
)
names(new_obs) <- c(par$output_obs_contamination, par$output_obs_clusters)
for (col in names(new_obs)) {
  anndataR::write_elem(
    par$output,
    file.path(obs_path, col),
    new_obs[[col]],
    overwrite = TRUE,
    compression = compression
  )
}

# New .obs columns are only read when listed in the 'column-order' attribute,
# written the same way as anndataR does
column_order <- union(
  as.character(rhdf5::h5readAttributes(par$output, obs_path)[["column-order"]]),
  names(new_obs)
)
output_h5 <- rhdf5::H5Fopen(par$output)
obs_group <- rhdf5::H5Gopen(output_h5, obs_path)
if (rhdf5::H5Aexists(obs_group, "column-order")) {
  rhdf5::H5Adelete(obs_group, "column-order")
}
rhdf5::h5writeAttribute(
  column_order, obs_group, "column-order",
  encoding = "UTF-8", variableLengthString = TRUE
)
rhdf5::H5Gclose(obs_group)
rhdf5::H5Fclose(output_h5)
