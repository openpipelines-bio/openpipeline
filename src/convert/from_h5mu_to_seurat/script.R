library(MuDataSeurat)
library(hdf5r)

## VIASH START
par <- list(
  input = "resources_test/10x_5k_anticmv/5k_human_antiCMV_T_TBNK_connect_mms.h5mu",
  output = "output.rds"
)
## VIASH END


tempfile <- tempfile(fileext=".h5mu")
file.copy(par$input, tempfile)

delete_modality <- function(open_h5, modality_path) {
  open_h5$link_delete(modality_path)
  mod_name <- sub("/mod/", "", modality_path)
  if ("mod-order" %in% names(h5attributes(open_h5[["mod"]]))) {
    current_attributes <- h5attributes(open_h5[["mod"]])$`mod-order`
    current_attributes <- current_attributes[current_attributes != mod_name] 
    h5attr(open_h5[["mod"]], "mod-order") <- current_attributes
  }
  for (obj_prefix in c("obsm/", "varm/", "varmap/", "obsmap/")) {
    obj_path = paste0(obj_prefix, mod_name)
    if (existsGroup(open_h5, obj_path)) {
      open_h5$link_delete(obj_path) 
    }
  }
  
} 

open_file <- H5File$new(tempfile, mode="r+")
modalities <- list.groups(open_file[["mod"]], full.names =  TRUE, recursive = FALSE)
to_delete = c()
for (modality_path in modalities) {
  dataset <- open_file[[modality_path]][["X"]]
  dataset_names <- names(dataset)
  if ("data" %in% dataset_names  && "indices" %in% dataset_names && "indptr" %in% dataset_names) {
    indexpointers <- dataset[["indptr"]]$read()
    indices <- dataset[["indices"]]$read()
    rowwise <- FALSE
    if ("encoding-type" %in% h5attr_names(dataset)) {
      rowwise <- h5attr(dataset, "encoding-type") == "csr_matrix"
    }
    if ("shape" %in% h5attr_names(dataset)) {
      X_dims <- h5attr(dataset, "shape")
    } else {
      X_dims <- c(length(indexpointers) - 1, max(indices) + 1)
      if (rowwise) {
        X_dims <- rev(X_dims)
      }
    }
    if (X_dims[2] < 1) {
      delete_modality(open_file, modality_path)
    }
  } else if (dataset$dims[1] < 1){
    delete_modality(open_file, modality_path)
  }
}

open_file$close_all()

cat("Reading input file\n")
obj <- ReadH5MU(tempfile)

cat("Writing output file\n")
saveRDS(obj, file = par$output, compress = TRUE)