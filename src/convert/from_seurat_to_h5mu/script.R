library(anndataR)
requireNamespace("reticulate", quietly = TRUE)
md <- reticulate::import("mudata")

### VIASH START
par <- list(
  input = "r./pbmc_1k_protein_v3_mms.rds",
  output = "./pbmc_1k_protein_v3_mms_converted.h5mu",
  assay = "RNA",
  modality = "rna",
  x_mapping = NULL,
  layers_mapping = TRUE,
  obs_mapping = TRUE,
  var_mapping = TRUE,
  uns_mapping = TRUE,
  obsm_mapping = TRUE,
  varm_mapping = TRUE,
  obsp_mapping = TRUE,
  varp_mapping = TRUE
)
### VIASH END

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
  assay_name = par$assay
)

# Create MuData object
mdata_obj <- md$MuData(setNames(list(h5ad_obj), par$modality))

# Save the MuData object as H5MU file
mdata_obj$write_h5mu(par$output)
