library(decontX)
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
    delta = c(10,10),
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

read_counts <- function(path, modality, layer) {
  mdata <- mudata$read_h5mu(path)
  adata <- mdata$mod[[modality]]
  list(
    mudata = mdata,
    adata = adata,
    # AnnData layers are backed by row-major (CSR) scipy sparse matrices,
    # which reticulate/anndata convert to R's dgRMatrix. decontX's internal
    # coercion to dgCMatrix does not support that class on recent Matrix
    # versions, so convert explicitly here via the generic sparse virtual
    # class.
    matrix = methods::as(Matrix::t(get_layer(adata, layer)), "CsparseMatrix")
  )
}

# Read input data
cat("Reading input file\n")
input <- read_counts(par$input, par$modality, par$input_layer)

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
  background <- read_counts(par$background, par$modality, background_layer)
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
# Store results back into the input AnnData object
input$adata$layers[[par$output_layer]] <- Matrix::t(result$decontXcounts)
input$adata$obs[[par$output_obs_contamination]] <- result$contamination
input$adata$obs[[par$output_obs_clusters]] <- result$z
# Writing output H5MU
input$mudata$write(par$output, compression = par$output_compression)