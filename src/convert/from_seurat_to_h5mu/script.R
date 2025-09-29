library(anndataR)
requireNamespace("reticulate", quietly = TRUE)
md <- reticulate::import("mudata")

### VIASH START
par <- list(
  input = "./pbmc_1k_protein_v3_mms.rds",
  output = "./pbmc_1k_protein_v3_mms_converted.h5mu",
  assay = "RNA",
  modality = "rna"
)
### VIASH END

seurat_obj <- readRDS(par$input)

h5ad_obj <- as_AnnData(
  seurat_obj,
  assay_name = par$assay
)

# Create MuData object
mdata_obj <- md$MuData(setNames(list(h5ad_obj), par$modality))

# Save the MuData object as H5MU file
mdata_obj$write_h5mu(par$output)
