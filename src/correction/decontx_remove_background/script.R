library(decontX)

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
# anndataR reads and writes .h5ad files but not .h5mu files, so hdf5r moves
# elements between both: the modality is extracted to a temporary .h5ad file
# for reading, and the results are written to a temporary .h5ad file with
# anndataR and copied into the output .h5mu.
# Files opened with hdf5r are closed with close_all(): objects left open keep
# the file locked for anndataR, which uses a different HDF5 library (rhdf5).

modality_exists <- function(h5file, modality_name) {
  h5file$exists("mod") && h5file[["mod"]]$exists(modality_name)
}

h5mu_to_h5ad <- function(h5mu_path, modality_name) {
  tmp_path <- tempfile(fileext = ".h5ad")
  src <- hdf5r::H5File$new(h5mu_path, mode = "r")
  on.exit(src$close_all())
  if (!modality_exists(src, modality_name)) {
    stop(
      "Modality '", modality_name, "' is not available in '", h5mu_path, "'."
    )
  }
  dest <- hdf5r::H5File$new(tmp_path, mode = "w")
  on.exit(dest$close_all(), add = TRUE)
  mod <- src[[paste0("mod/", modality_name)]]

  for (child in mod$names) {
    dest$obj_copy_from(mod, child, child)
  }
  # The modality group's attributes (encoding-type: anndata, ...) become the
  # root attributes of the .h5ad file
  for (name in hdf5r::h5attr_names(mod)) {
    value <- hdf5r::h5attr(mod, name)
    space <- if (length(value) == 1) {
      hdf5r::H5S$new("scalar")
    } else {
      hdf5r::guess_space(value)
    }
    dest$create_attr(name, value, space = space)
  }

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

# Replace (or add) element 'name' of 'dest_group' with 'src_name' of 'src'
h5_replace <- function(dest_group, name, src, src_name) {
  if (dest_group$exists(name)) dest_group$link_delete(name)
  dest_group$obj_copy_from(src, src_name, name)
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
# Write only the results as an AnnData object to a temporary .h5ad file
obs <- data.frame(row.names = colnames(input$matrix))
obs[[par$output_obs_contamination]] <- as.numeric(result$contamination)
obs[[par$output_obs_clusters]] <- factor(result$z)
layers <- list()
layers[[par$output_layer]] <- methods::as(
  Matrix::t(result$decontXcounts), "RsparseMatrix"
)
results <- anndataR::AnnData(
  obs = obs,
  var = data.frame(row.names = rownames(input$matrix)),
  layers = layers
)
results_path <- tempfile(fileext = ".h5ad")
compression <- if (is.null(par$output_compression)) {
  "none"
} else {
  par$output_compression
}
anndataR::write_h5ad(results, results_path, compression = compression)

# Copy the input and add the results to it in place, leaving every other
# modality and field untouched.
if (!file.copy(par$input, par$output, overwrite = TRUE)) {
  stop("Could not copy '", par$input, "' to '", par$output, "'.")
}
output <- hdf5r::H5File$new(par$output, mode = "r+")
results_h5 <- hdf5r::H5File$new(results_path, mode = "r")
mod <- output[[paste0("mod/", par$modality)]]

if (!mod$exists("layers")) mod$create_group("layers")
h5_replace(
  mod[["layers"]], par$output_layer,
  results_h5, paste0("layers/", par$output_layer)
)

obs_group <- mod[["obs"]]
new_obs_columns <- c(par$output_obs_contamination, par$output_obs_clusters)
for (col in new_obs_columns) {
  h5_replace(obs_group, col, results_h5, paste0("obs/", col))
}

# New .obs columns are only read when listed in the 'column-order' attribute
column_order <- union(
  as.character(hdf5r::h5attr(obs_group, "column-order")), new_obs_columns
)
obs_group$attr_delete("column-order")
obs_group$create_attr("column-order", column_order)

results_h5$close_all()
output$close_all()
unlink(results_path)
