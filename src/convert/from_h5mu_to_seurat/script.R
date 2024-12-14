library(MuDataSeurat)
library(hdf5r)

## VIASH START
par <- list(
  input = paste0(
    "resources_test/10x_5k_anticmv/",
    "5k_human_antiCMV_T_TBNK_connect_mms.h5mu"
  ),
  output = "output.rds"
)
## VIASH END


temp_h5mu <- tempfile(fileext = ".h5mu")
file.copy(par$input, temp_h5mu)

delete_modality <- function(open_h5, modality_path) {
  open_h5$link_delete(modality_path)
  mod_name <- sub("/mod/", "", modality_path)
  if ("mod-order" %in% names(hdf5r::h5attributes(open_h5[["mod"]]))) {
    current_attributes <- hdf5r::h5attributes(open_h5[["mod"]])$`mod-order`
    current_attributes <- current_attributes[current_attributes != mod_name]
    hdf5r::h5attr(open_h5[["mod"]], "mod-order") <- current_attributes
  }
  for (obj_prefix in c("obsm/", "varm/", "varmap/", "obsmap/")) {
    obj_path <- paste0(obj_prefix, mod_name)
    if (hdf5r::existsGroup(open_h5, obj_path)) {
      open_h5$link_delete(obj_path)
    }
  }
}

open_file <- H5File$new(temp_h5mu, mode = "r+")
modalities <- list.groups(open_file[["mod"]],
  full.names = TRUE, recursive = FALSE
)
to_delete <- c()

determine_matrix_dims <- function(dataset, indexpointers, indices, rowwise) {
  if ("shape" %in% hdf5r::h5attr_names(dataset)) {
    x_dims <- hdf5r::h5attr(dataset, "shape")
  } else {
    x_dims <- c(length(indexpointers) - 1, max(indices) + 1)
    if (rowwise) {
      x_dims <- rev(x_dims)
    }
  }
  return(x_dims)
}
for (modality_path in modalities) {
  dataset <- open_file[[modality_path]][["X"]]
  dataset_names <- names(dataset)
  if (
    "data" %in% dataset_names &&
      "indices" %in% dataset_names &&
      "indptr" %in% dataset_names
  ) {
    indexpointers <- dataset[["indptr"]]$read()
    indices <- dataset[["indices"]]$read()
    rowwise <- FALSE
    if ("encoding-type" %in% h5attr_names(dataset)) {
      rowwise <- h5attr(dataset, "encoding-type") == "csr_matrix"
    }
    x_dims <- determine_matrix_dims(dataset, indexpointers, indices, rowwise)
    if (x_dims[2] < 1) {
      delete_modality(open_file, modality_path)
    }
  } else if (dataset$dims[1] < 1) {
    delete_modality(open_file, modality_path)
  }
}

open_file$close_all()

cat("Reading input file\n")
obj <- ReadH5MU(temp_h5mu)

cat("Writing output file\n")
saveRDS(obj, file = par$output, compress = TRUE)
