library(anndataR)

rhdf5::h5disableFileLocking()

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

# Reading a single modality from a MuData file via `group` requires
# anndataR >= 1.4.0
stopifnot(packageVersion("anndataR") >= "1.4.0")

# If the input is a MuData file, read a single modality as input instead
if (is_mudata_file(par$input)) {
  if (is.null(par$modality) || par$modality == "") {
    stop("'modality' argument must be set if the input is a MuData file.")
  }
  group <- paste0("mod/", par$modality)
} else {
  group <- "/"
}

seurat_obj <- read_h5ad(
  par$input,
  mode = "r",
  as = "Seurat",
  group = group,
  assay_name = par$assay,
  x_mapping = par$x_mapping
)

saveRDS(seurat_obj, file = par$output)
