### VIASH START
par <- list(
  input = "test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5ad",
  output = "output.rds"
)
### VIASH END

library(SeuratDisk)

print(par$input)

# SeuratDisk needs to convert to h5seurat first
tempfile <- paste0(gsub("\\..*$", "", par$output), ".tmp.h5seurat")
if (file.exists(tempfile)) file.remove(tempfile)

# remove temporary h5seurat file after conversion
on.exit(file.remove(tempfile))

# connect to input file
# ignore the warning on "unknown file type: h5ad"
# cat("connect to input file\n")
hfile <- Connect(par$input, type = "h5ad", force = TRUE, mode = "r")

# convert h5ad to h5seurat
# cat("convert h5ad to h5seurat\n")
converted <- Convert(hfile, dest = tempfile, overwrite = TRUE)

# load h5seurat
# cat("load h5seurat\n")
seuratObject <- LoadH5Seurat(converted)

# save as rds
saveRDS(seuratObject, file = par$output, compress = TRUE)
