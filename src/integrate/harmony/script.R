library(reticulate)
library(harmony)

mudata <- reticulate::import("mudata")

### VIASH START
par <- list(
  input = "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
  theta = 0,
  obs_name = "batch"
)
### VIASH END

data <- mudata$read_h5mu(par$input)
rna_data <- data$mod[["rna"]]
pca_embedding <- rna_data$obsm[["X_pca"]]
meta_data <- rna_data$obs

## Run Harmonydata_mat, meta_data, vars_use
harmony_embedding <- HarmonyMatrix(
  data_mat = pca_embedding,
  meta_data = meta_data,
  vars_use = par$obs_name,
  do_pca = FALSE,
  theta = par$theta
)

## Add Harmony embeddings to Anndata
rna_data$obsm[["X_harmony"]] <- harmony_embedding

## Save as H5AD
data$write(par$output)