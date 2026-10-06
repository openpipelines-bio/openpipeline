library(anndataR)
library(Azimuth)
library(Seurat)
requireNamespace("hdf5r", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
md <- reticulate::import("mudata")

rhdf5::h5disableFileLocking()

### VIASH START
par <- list(
  input = paste0(
    "resources_test/pbmc_1k_protein_v3/",
    "pbmc_1k_protein_v3_mms.rds"
  ),
  modality = "rna",
  assay = "RNA",
  reference = "pbmcref",
  annotation_levels = NULL,
  do_adt = FALSE,
  umap_name = "ref.umap",
  k_weight = 50,
  n_trees = 20,
  mapping_score_k = 100,
  output = "output.rds",
  output_compression = NULL
)
### VIASH END

# Read input data
cat("Reading input file\n")
query <- readRDS(par$input)


cat("Running Azimuth reference mapping against reference:", par$reference, "\n")
obj <- RunAzimuth(
  query = query,
  reference = par$reference,
  annotation.levels = par$annotation_levels,
  umap.name = par$umap_name,
  do.adt = par$do_adt,
  assay = par$assay,
  k.weight = par$k_weight,
  n.trees = par$n_trees,
  mapping.score.k = par$mapping_score_k,
  verbose = TRUE
)

cat("Writing output data\n")

# Ouput can be in Seurat format
cat("Writing output file\n")
saveRDS(obj, file = par$output, compress = par$output_compression)
