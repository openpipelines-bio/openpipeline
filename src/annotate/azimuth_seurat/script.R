library(Azimuth)
library(Seurat)

### VIASH START
par <- list(
  input = paste0(
    "resources_test/pbmc_1k_protein_v3/",
    "pbmc_1k_protein_v3_mms.rds"
  ),
  assay = "RNA",
  input_layer = "X",
  input_var_gene_names = "gene_symbol",
  reference = "pbmcref",
  annotation_levels = NULL,
  do_adt = FALSE,
  umap_name = "ref.umap",
  k_weight = 50,
  n_trees = 20,
  mapping_score_k = 100,
  output = "output.rds",
  output_compression = "gzip"
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

# Run Azimuth on a minimal query object containing just the raw counts,
# and transfer the results back onto the original object afterwards.
cat("Creating query object from layer '", par$input_layer, "' of assay '",
  par$assay, "'\n",
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
query <- CreateSeuratObject(counts = counts, assay = "RNA")
rm(counts)

cat("Running Azimuth reference mapping against reference:", par$reference, "\n")
query <- RunAzimuth(
  query = query,
  reference = par$reference,
  annotation.levels = par$annotation_levels,
  umap.name = par$umap_name,
  do.adt = par$do_adt,
  assay = "RNA",
  k.weight = par$k_weight,
  n.trees = par$n_trees,
  mapping.score.k = par$mapping_score_k,
  verbose = TRUE
)

cells <- Cells(obj)
if (!setequal(Cells(query), cells)) {
  stop("Cells of the mapped query do not match the cells of the input object.")
}

warn_overwrite <- function(name, existing, what) {
  if (name %in% existing) {
    warning(what, " '", name, "' already exists in the input object ",
      "and will be overwritten.",
      call. = FALSE
    )
  }
}

cat("Adding predicted labels and mapping scores to the cell metadata\n")
result_cols <- grep(
  "^predicted\\.|^mapping\\.score$", colnames(query[[]]),
  value = TRUE
)
for (col in result_cols) {
  warn_overwrite(col, colnames(obj[[]]), "Metadata column")
}
obj <- AddMetaData(obj, metadata = query[[]][cells, result_cols, drop = FALSE])

# Cell embeddings are not guaranteed to be in the same order as the cells of
# the object after IntegrateEmbeddings(), so realign them to the input object.
cat("Adding reference-projected reductions\n")
for (reduction_name in c(par$umap_name, "integrated_dr")) {
  warn_overwrite(reduction_name, Reductions(obj), "Reduction")
  reduction <- query[[reduction_name]]
  reduction <- CreateDimReducObject(
    embeddings = Embeddings(reduction)[cells, , drop = FALSE],
    loadings = Loadings(reduction),
    key = paste0(gsub("[^[:alnum:]]", "", reduction_name), "_"),
    assay = par$assay,
    misc = Misc(reduction)
  )
  obj[[reduction_name]] <- reduction
}

# TransferData() stores the per-class prediction scores of each annotation
# level (and, if requested, the imputed ADT expression) as separate assays.
cat("Adding prediction score assays\n")
for (assay_name in setdiff(Assays(query), "RNA")) {
  warn_overwrite(assay_name, Assays(obj), "Assay")
  data <- LayerData(query, assay = assay_name, layer = "data")
  obj[[assay_name]] <- CreateAssayObject(data = data[, cells, drop = FALSE])
}

cat("Writing output file\n")
compress <- if (par$output_compression == "none") {
  FALSE
} else {
  par$output_compression
}
saveRDS(obj, file = par$output, compress = compress)
