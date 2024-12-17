library(reticulate)
library(harmony)
library(anndata)

mudata <- reticulate::import("mudata")

### VIASH START
par <- list(
  input = "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
  output = "foo.h5mu",
  modality = "rna",
  obsm_input = "X_pca",
  obsm_output = "X_pca_int",
  theta = 2,
  obs_covariates = c("leiden")
)
### VIASH END

data <- mudata$read_h5mu(par$input)
rna_data <- data$mod[[par$modality]]

## Run Harmony
pca_in <- rna_data$obsm[[par$obsm_input]]
meta_data <- rna_data$obs
harmony_embedding <- HarmonyMatrix(
  data_mat = pca_in,
  meta_data = meta_data,
  vars_use = par$obs_covariates,
  theta = par$theta,
  do_pca = FALSE
)

## Add Harmony embeddings to Anndata
rna_data$obsm[[par$obsm_output]] <- harmony_embedding

## Save as h5mu
data$write(par$output, compression = par$output_compression)
