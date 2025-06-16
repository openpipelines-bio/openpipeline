library(anndataR)

### VIASH START
par <- list(
  input = "data_split/split_data.h5ad",
  output = "data_split/split_data_anndatar.rds",
  assay = "RNA"
)
### VIASH END


seurat_obj <- read_h5ad(
  par$input,
  mode = "r",
  as = "Seurat",
  assay_name = par$assay
)

saveRDS(seurat_obj, file = par$output)
