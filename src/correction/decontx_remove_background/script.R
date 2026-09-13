library(decontX)
library(anndataR)
# Loaded (but not called directly) so reticulate's R<->Python sparse-matrix
# converters are registered for the write-back step below, which still goes
# through Python's mudata via reticulate (anndataR has no multi-modal .h5mu
# writer, so it's only used for reading here).
requireNamespace("anndata", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
mudata <- reticulate::import("mudata")

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

get_layer <- function(adata, layer) {
  # find data
  data <-
    if (is.null(layer)) {
      adata$X
    } else {
      adata$layers[[layer]]
    }

  # check if data is available
  if (is.null(data)) {
    if (is.null(layer)) {
      stop("No layer specified and no .X slot available in the AnnData object.")
    } else {
      stop(
        "Requested layer '",
        layer,
        "' is not available in the AnnData object. Available layers: ",
        paste(names(adata$layers), collapse = ", ")
      )
    }
  }

  # Set matrix dimnames
  dimnames(data) <- list(adata$obs_names, adata$var_names)

  # return output
  data
}

# Read a single modality's counts + obs with anndataR, given an
# already-loaded (reticulate) Python AnnData object for that modality.
# anndataR only reads standalone .h5ad files, so the modality is dumped there
# via its own (cheap, in-memory) write_h5ad() rather than re-reading the
# source .h5mu file from disk a second time.
read_counts <- function(py_mod, layer) {
  h5ad_path <- tempfile(fileext = ".h5ad")
  py_mod$write_h5ad(h5ad_path)
  adata <- anndataR::read_h5ad(h5ad_path, as = "InMemoryAnnData")
  list(
    adata = adata,
    # decontX's internal coercion to dgCMatrix does not support every sparse
    # subclass on recent Matrix versions, so convert explicitly here via the
    # generic sparse virtual class.
    matrix = methods::as(Matrix::t(get_layer(adata, layer)), "CsparseMatrix")
  )
}

# Read input data
cat("Reading input file\n")
input_mdata <- mudata$read_h5mu(par$input)
input <- read_counts(input_mdata$mod[[par$modality]], par$input_layer)

# Read background data, if provided
background_matrix <- NULL
background_batch <- NULL
if (!is.null(par$background)) {
  cat("Reading background file\n")
  background_layer <- if (is.null(par$background_layer)) {
    par$input_layer
  } else {
    par$background_layer
  }
  background_mdata <- mudata$read_h5mu(par$background)
  background <- read_counts(
    background_mdata$mod[[par$modality]], background_layer
  )
  background_matrix <- background$matrix
  if (!is.null(par$background_obs_batch)) {
    background_batch <- background$adata$obs[[par$background_obs_batch]]
  }
}

# Get optional cluster / batch labels from .obs
input_clusters <- if (is.null(par$input_obs_clusters)) {
  NULL
} else {
  input$adata$obs[[par$input_obs_clusters]]
}
input_batch <- if (is.null(par$input_obs_batch)) {
  NULL
} else {
  input$adata$obs[[par$input_obs_batch]]
}

cat("Estimating and removing contamination with DecontX\n")
result <- decontX(
  x = input$matrix,
  z = input_clusters,
  batch = input_batch,
  background = background_matrix,
  bgBatch = background_batch,
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
# Store results back into the already-loaded input MuData object (no second
# read of par$input) and write the full MuData out, preserving every other
# modality and field untouched.
output_mod <- input_mdata$mod[[par$modality]]
output_mod$layers[[par$output_layer]] <- Matrix::t(result$decontXcounts)
output_mod$obs[[par$output_obs_contamination]] <- result$contamination
output_mod$obs[[par$output_obs_clusters]] <- result$z
input_mdata$write(par$output, compression = par$output_compression)
