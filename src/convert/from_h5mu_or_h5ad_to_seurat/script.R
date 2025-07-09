library(anndataR)
library(hdf5r)


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
  tmp_path <- tempfile(fileext = ".h5ad")
  mod_location <- paste("mod", modality_name, sep = "/")
  h5src <- hdf5r::H5File$new(h5mu_path, "r")
  h5dest <- hdf5r::H5File$new(tmp_path, "w")
  # Copy over the child objects and the child attributes from root
  # Root cannot be copied directly because it always exists and
  # copying does not allow overwriting.
  children <- hdf5r::list.objects(h5src,
    path = mod_location,
    full.names = FALSE, recursive = FALSE
  )
  for (child in children) {
    h5dest$obj_copy_from(
      h5src, paste(mod_location, child, sep = "/"),
      paste0("/", child)
    )
  }
  # Also copy the root attributes
  root_attrs <- hdf5r::h5attr_names(x = h5src)
  for (attr in root_attrs) {
    h5a <- h5src$attr_open(attr_name = attr)
    robj <- h5a$read()
    h5dest$create_attr_by_name(
      attr_name = attr,
      obj_name = ".",
      robj = robj,
      space = h5a$get_space(),
      dtype = h5a$get_type()
    )
  }
  h5src$close()
  h5dest$close()

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
