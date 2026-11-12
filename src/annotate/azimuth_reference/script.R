library(Azimuth)
library(Seurat)

### VIASH START
par <- list(
  input = "resources_test/annotation_test_data/TS_Blood_filtered.rds",
  assay = "RNA",
  input_layer = "counts",
  input_var_gene_names = NULL,
  annotation_levels = c("cell_ontology_class", "cell_type"),
  reference_version = "0.0.0",
  n_variable_features = 3000,
  n_dims = 50,
  k_param = 31,
  seed = 42,
  output = "azimuth_reference"
)
### VIASH END

cat("Reading input file\n")
obj <- readRDS(par$input)
if (!inherits(obj, "Seurat")) {
  stop("Input file should contain a Seurat object, found: ", class(obj)[[1]])
}
if (!par$assay %in% Assays(obj)) {
  stop(
    "Assay '", par$assay, "' not found in input object. Available assays: ",
    paste(Assays(obj), collapse = ", ")
  )
}
if (!par$input_layer %in% Layers(obj[[par$assay]])) {
  stop(
    "Layer '", par$input_layer, "' not found in assay '", par$assay,
    "'. Available layers: ", paste(Layers(obj[[par$assay]]), collapse = ", ")
  )
}

missing_levels <- setdiff(par$annotation_levels, colnames(obj[[]]))
if (length(missing_levels) > 0) {
  stop(
    "Annotation level(s) not found in the cell metadata of the input object: ",
    paste(missing_levels, collapse = ", "), ". Available columns: ",
    paste(colnames(obj[[]]), collapse = ", ")
  )
}
if (anyDuplicated(par$annotation_levels)) {
  stop("--annotation_levels should be unique.")
}
# AzimuthReference() stores the "ori.index" column next to the annotation
# levels, and RunAzimuth() skips levels starting with "ori", "nCount" or
# "nFeature" when transferring all levels.
reserved <- grepl("^(ori|nCount|nFeature)", par$annotation_levels)
if (any(reserved)) {
  stop(
    "Annotation levels should not start with 'ori', 'nCount' or 'nFeature', ",
    "found: ", paste(par$annotation_levels[reserved], collapse = ", ")
  )
}

if (par$n_dims > ncol(obj)) {
  stop(
    "--n_dims (", par$n_dims, ") should not be larger than the number of ",
    "cells in the input object (", ncol(obj), ")."
  )
}

# Build the reference from a minimal object containing just the raw counts
# and the annotation levels.
cat("Creating reference object from layer '", par$input_layer,
  "' of assay '", par$assay, "'\n",
  sep = ""
)
counts <- LayerData(obj, assay = par$assay, layer = par$input_layer)
if (!is.null(par$input_var_gene_names)) {
  feature_meta <- obj[[par$assay]][[]]
  if (!par$input_var_gene_names %in% colnames(feature_meta)) {
    stop(
      "Column '", par$input_var_gene_names, "' not found in the feature ",
      "metadata of assay '", par$assay, "'. Available columns: ",
      paste(colnames(feature_meta), collapse = ", ")
    )
  }
  gene_names <- as.character(
    feature_meta[rownames(counts), par$input_var_gene_names]
  )
  # Fall back to the original feature name where no gene name is available
  missing <- is.na(gene_names) | gene_names == ""
  gene_names[missing] <- rownames(counts)[missing]
  rownames(counts) <- make.unique(gene_names)
}

metadata <- obj[[]][colnames(counts), par$annotation_levels, drop = FALSE]
for (level in par$annotation_levels) {
  labels <- metadata[[level]]
  if (any(is.na(labels))) {
    stop(
      "Annotation level '", level, "' contains ", sum(is.na(labels)),
      " cell(s) without a label. Remove these cells from the input object, ",
      "or assign them a label (e.g. 'Unknown')."
    )
  }
  # Unused factor levels would end up in the reference color map, which
  # AzimuthReference() requires to match the labels exactly.
  metadata[[level]] <- if (is.factor(labels)) {
    droplevels(labels)
  } else {
    factor(labels, levels = sort(unique(as.character(labels))))
  }
}
ref <- CreateSeuratObject(counts = counts, assay = "RNA", meta.data = metadata)
rm(obj, counts, metadata)

cat("Normalizing counts using SCTransform\n")
ref <- SCTransform(
  ref,
  assay = "RNA",
  new.assay.name = "SCT",
  variable.features.n = par$n_variable_features,
  seed.use = par$seed,
  verbose = TRUE
)

# The supervised PCA is guided by the neighbor graph of a regular PCA
cat("Computing supervised PCA\n")
dims <- seq_len(par$n_dims)
ref <- RunPCA(ref, assay = "SCT", npcs = par$n_dims, verbose = FALSE)
ref <- FindNeighbors(
  ref,
  reduction = "pca",
  dims = dims,
  graph.name = c("pca_nn", "pca_snn"),
  verbose = FALSE
)
ref <- RunSPCA(
  ref,
  assay = "SCT",
  graph = "pca_snn",
  npcs = par$n_dims,
  reduction.name = "spca",
  verbose = FALSE
)

# RunAzimuth() projects the query onto the reference UMAP, which requires
# the UMAP model to be stored.
cat("Computing UMAP model\n")
ref <- RunUMAP(
  ref,
  reduction = "spca",
  dims = dims,
  reduction.name = "umap",
  return.model = TRUE,
  seed.use = par$seed,
  verbose = FALSE
)

cat("Creating Azimuth reference\n")
ref <- AzimuthReference(
  ref,
  refUMAP = "umap",
  refDR = "spca",
  refAssay = "SCT",
  dims = dims,
  k.param = par$k_param,
  plotref = "umap",
  metadata = par$annotation_levels,
  reference.version = par$reference_version,
  verbose = TRUE
)

cat("Writing reference to", par$output, "\n")
dir.create(par$output, recursive = TRUE, showWarnings = FALSE)
SaveAzimuthReference(
  ref,
  folder = paste0(normalizePath(par$output), "/")
)
