library(anndataR)
library(Azimuth)
library(Seurat)
requireNamespace("hdf5r", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
md <- reticulate::import("mudata")

rhdf5::h5disableFileLocking()

### VIASH START
par <- list(
  input = "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
  modality = "rna",
  assay = "RNA",
  reference = "pbmcref",
  annotation_levels = NULL,
  do_adt = FALSE,
  umap_name = "ref.umap",
  k_weight = 50,
  n_trees = 20,
  mapping_score_k = 100,
  output = "output.h5mu",
  output_compression = NULL
)
### VIASH END

# Copy a single modality out of a MuData (h5mu) file into a standalone
# h5ad file, so it can be read as a plain AnnData/Seurat object. Adapted
# from convert/from_h5mu_or_h5ad_to_seurat/script.R.
h5mu_to_h5ad <- function(h5mu_path, modality_name) {
  tmp_path <- tempfile(fileext = ".h5ad")
  mod_location <- paste("mod", modality_name, sep = "/")
  h5src <- hdf5r::H5File$new(h5mu_path, "r")
  h5dest <- hdf5r::H5File$new(tmp_path, "w")
  children <- hdf5r::list.objects(h5src,
    path = mod_location,
    full.names = FALSE, recursive = FALSE
  )
  for (child in children) {
    h5dest$obj_copy_from(
      h5src, paste(mod_location, child, sep = "/"),
      paste0("/", child)
    )
  }
  root_attrs <- hdf5r::h5attr_names(x = h5src)
  for (attr in root_attrs) {
    h5a <- h5src$attr_open(attr_name = attr)
    robj <- h5a$read()
    h5dest$create_attr_by_name(
      attr_name = attr,
      obj_name = ".",
      robj = robj,
      space = h5a$get_space(),
      dtype = h5a$get_type()
    )
  }
  h5src$close()
  h5dest$close()

  tmp_path
}

cat("Reading input file\n")
h5ad_path <- h5mu_to_h5ad(par$input, par$modality)
query <- read_h5ad(
  h5ad_path,
  mode = "r",
  as = "Seurat",
  assay_name = par$assay,
  x_mapping = "counts"
)

cat("Running Azimuth reference mapping against reference:", par$reference, "\n")
query <- RunAzimuth(
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

# RunAzimuth()'s integration/projection steps (e.g. IntegrateEmbeddings())
# do not guarantee that a reduction's cell.embeddings rows stay in the same
# order as the parent object's cells. anndataR strictly requires obsm row
# names to match .obs_names exactly (including order), so realign every
# reduction here before converting.
for (reduction_name in Reductions(query)) {
  embeddings <- Embeddings(query[[reduction_name]])
  query[[reduction_name]]@cell.embeddings <- embeddings[Cells(query), , drop = FALSE]
}

annotated_h5ad <- as_AnnData(
  query,
  assay_name = par$assay,
  output_class = c("ReticulateAnnData")
)

# Read the full original MuData object so that any other modalities
# (e.g. "prot") are preserved untouched in the output. RunAzimuth()
# converts query gene IDs to gene symbols AND subsets to the overlapping
# feature set, so the annotated modality's var_names are simultaneously
# renamed and filtered relative to the original. Mutating an existing
# MuData's modality in place in that situation isn't supported: MuData
# can't reconcile its cached global var index and raises
# "var_names seem to have been renamed and filtered at the same time",
# explicitly suggesting to instead build a fresh MuData from the
# modalities dict, which is what's done here.
input_mudata <- md$read_h5mu(par$input)
mod_names <- reticulate::iterate(input_mudata$mod$keys())
mods <- reticulate::dict()
for (mod_name in mod_names) {
  mods[[mod_name]] <- if (identical(mod_name, par$modality)) {
    annotated_h5ad
  } else {
    input_mudata$mod[[mod_name]]
  }
}
output_mudata <- md$MuData(mods)
output_mudata$write_h5mu(par$output, compression = par$output_compression)
