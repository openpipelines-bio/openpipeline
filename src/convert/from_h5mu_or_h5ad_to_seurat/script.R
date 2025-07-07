library(anndataR)
library(hdf5r)
library(reticulate)

mudata <- reticulate::import("mudata")
anndata <- reticulate::import("anndata")

### VIASH START
par <- list(
  input = "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
  output = "resources_test/annotation_test_data/TS_Blood_filtered.rds",
  assay = "RNA",
  modality = "rna"
)
### VIASH END

is_mudata_file <- function(file_path) {
  con <- file(file_path, "rb")
  tryCatch({
    # It is possible to create a MuData compatible h5 file without using
    # MuData, and those could not have "MuData" as the first bytes.
    # But that is really an edge case and this check should hold.
    header <- readBin(con, "raw", n = 6)
    identical(header, charToRaw("MuData"))
  }, finally = {
    close(con)
  })
}

h5mu_to_h5ad <- function(h5mu_path, modality_name) {
  # Conversion of h5mu file to h5ad file,
  # using python's mudata/anndata with reticulate
  tmp_path <- tempfile(fileext = ".h5ad")
  adata <- mudata$read_h5ad(par$input, par$modality)
  adata$write_h5ad(tmp_path, compression = "gzip")
  tmp_path
}

# If the input is a MuData file, use a single modality as input instead
if (is_mudata_file(par$input)) {
  if (is.null(par$modality) || par$modality == "") {
    stop("'modality' argument must be set if the input is a MuData file.")
  }
  h5ad_path <- h5mu_to_h5ad(par$input, par$modality)
} else {
  h5ad_path <- par$input
}

seurat_obj <- read_h5ad(
  h5ad_path,
  mode = "r",
  as = "Seurat",
  assay_name = par$assay
)

saveRDS(seurat_obj, file = par$output)
