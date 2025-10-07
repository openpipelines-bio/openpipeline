library(anndataR)
requireNamespace("reticulate", quietly = TRUE)
md <- reticulate::import("mudata")

### VIASH START
par <- list(
  input = "./pbmc_1k_protein_v3_mms.rds",
  output = "./pbmc_1k_protein_v3_mms_converted.h5mu",
  assay = "RNA",
  modality = "rna",
  x_mapping = NULL,
  layers_mapping = "True",
  obs_mapping = "True",
  var_mapping = "True",
  uns_mapping = "True",
  obsm_mapping = "True",
  varm_mapping = "True",
  obsp_mapping = "True",
  varp_mapping = "True"
  modality = "rna",
  x_mapping = NULL,
  layers_mapping = "True",
  obs_mapping = "True",
  var_mapping = "True",
  uns_mapping = "True",
  obsm_mapping = "True",
  varm_mapping = "True",
  obsp_mapping = "True",
  varp_mapping = "True"
)
### VIASH END

process_mapping_param <- function(param_key) {
  param_value <- par[[param_key]]
  if (is.null(param_value)) {
    NULL
  } else if (tolower(param_value) == "true") {
    TRUE
  } else if (tolower(param_value) == "false") {
    FALSE
  } else {
    # Try to parse as JSON for named character vector
    tryCatch(
      {
        return(jsonlite::fromJSON(param_value))
      },
      error = function(e) {
        cat("Could not parse json from argument", param_key, "\n")
        cat("Provided value:", param_value, "\n")
        stop(e)
      }
    )
  }
}

mapping_values <- c(
  "layers_mapping",
  "obs_mapping",
  "var_mapping",
  "uns_mapping",
  "obsm_mapping",
  "varm_mapping",
  "obsp_mapping",
  "varp_mapping"
)

for (mapping in mapping_values) {
  par[[mapping]] <- process_mapping_param(mapping)
}

process_mapping_param <- function(param_key) {
  param_value <- par[[param_key]]
  if (is.null(param_value)) {
    NULL
  } else if (tolower(param_value) == "true") {
    TRUE
  } else if (tolower(param_value) == "false") {
    FALSE
  } else {
    # Try to parse as JSON for named character vector
    tryCatch(
      {
        return(jsonlite::fromJSON(param_value))
      },
      error = function(e) {
        cat("Could not parse json from argument", param_key, "\n")
        cat("Provided value:", param_value, "\n")
        stop(e)
      }
    )
  }
}

mapping_values <- c(
  "layers_mapping",
  "obs_mapping",
  "var_mapping",
  "uns_mapping",
  "obsm_mapping",
  "varm_mapping",
  "obsp_mapping",
  "varp_mapping"
)

for (mapping in mapping_values) {
  par[[mapping]] <- process_mapping_param(mapping)
}

seurat_obj <- readRDS(par$input)
h5ad_obj <- as_AnnData(
  seurat_obj,
  layers_mapping = par$layers_mapping,
  obs_mapping = par$obs_mapping,
  var_mapping = par$var_mapping,
  uns_mapping = par$uns_mapping,
  obsm_mapping = par$obsm_mapping,
  varm_mapping = par$varm_mapping,
  obsp_mapping = par$obsp_mapping,
  varp_mapping = par$varp_mapping,
  x_mapping = par$x_mapping,
  assay_name = par$assay,
  output_class = "ReticulateAnnData"
)

# Create MuData object
mods <- reticulate::dict()
mods[[par$modality]] <- h5ad_obj
mdata_obj <- md$MuData(mods)

# Save the MuData object as H5MU file
mdata_obj$write_h5mu(par$output, compression = par$output_compression)
