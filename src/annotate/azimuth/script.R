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
  input_reference_gene_overlap = 100,
  reference = "pbmcref",
  custom_reference = NULL,
  annotation_levels = NULL,
  do_adt = FALSE,
  k_weight = 50,
  n_trees = 20,
  mapping_score_k = 100,
  output = "output.rds",
  output_compression = "gzip",
  output_obs_predictions = NULL,
  output_obs_probability = NULL,
  output_obs_mapping_score = "azimuth_mapping_score",
  output_obsm_umap = "X_azimuth_umap",
  output_obsm_integrated = "X_integrated_azimuth"
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

if (par$reference == "custom") {
  if (is.null(par$custom_reference)) {
    stop("--custom_reference must be provided when --reference is 'custom'.")
  }
  if (!dir.exists(par$custom_reference)) {
    stop("--custom_reference should be an existing directory, found: ",
      par$custom_reference,
      call. = FALSE
    )
  }
  missing_files <- setdiff(
    c("ref.Rds", "idx.annoy"),
    list.files(par$custom_reference)
  )
  if (length(missing_files) > 0) {
    stop(
      "--custom_reference is missing the required file(s): ",
      paste(missing_files, collapse = ", "),
      call. = FALSE
    )
  }
  reference <- par$custom_reference
} else {
  if (!is.null(par$custom_reference)) {
    stop("--custom_reference can only be used with --reference 'custom'.")
  }
  reference <- par$reference
}

# The output names have to line up with the annotation levels, which are
# only known when --annotation_levels is provided.
for (arg in c("output_obs_predictions", "output_obs_probability")) {
  if (is.null(par[[arg]])) {
    next
  }
  if (is.null(par$annotation_levels)) {
    stop("--", arg, " can only be used together with --annotation_levels.")
  }
  if (length(par[[arg]]) != length(par$annotation_levels)) {
    stop(
      "--", arg, " should have the same length as --annotation_levels (",
      length(par$annotation_levels), "), found: ", length(par[[arg]]), "."
    )
  }
}
if (par$output_obsm_umap == par$output_obsm_integrated) {
  stop("--output_obsm_umap and --output_obsm_integrated should be different.")
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

# Check overlap of genes between the query and the reference
cat("Checking overlap of genes between input and reference\n")
# Load the reference the same way RunAzimuth() does: custom references from
# a directory, curated references via SeuratData (installed on first use).
if (par$reference == "custom") {
  ref_map <- LoadReference(reference)$map
} else {
  if (!reference %in% SeuratData::InstalledData()$Dataset) {
    SeuratData::InstallData(reference)
  }
  ref_map <- SeuratData::LoadData(reference, type = "azimuth")$map
}
ref_genes <- rownames(ref_map)
rm(ref_map)
# RunAzimuth() converts the query gene names (e.g. Ensembl IDs) to the
# reference gene symbols, so do the same before checking the overlap.
# RunAzimuth() leaves the already converted names untouched.
query <- ConvertGeneNames(
  object = query,
  reference.names = ref_genes,
  homolog.table = "https://seurat.nygenome.org/azimuth/references/homologs.rds"
)
common_genes <- intersect(ref_genes, rownames(query))
if (length(common_genes) < par$input_reference_gene_overlap) {
  stop(
    "The intersection of genes between the query and reference is too small.\n",
    "Expected overlap: ", par$input_reference_gene_overlap, "\n",
    "Detected overlap: ", length(common_genes)
  )
}
rm(ref_genes, common_genes)

cat("Running Azimuth reference mapping against reference:", reference, "\n")
query <- RunAzimuth(
  query = query,
  reference = reference,
  annotation.levels = par$annotation_levels,
  umap.name = "ref.umap",
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

# RunAzimuth() stores its results as `predicted.<level>`,
# `predicted.<level>.score` and `mapping.score`. When no annotation levels
# were requested, the transferred levels are derived from these columns.
annotation_levels <- par$annotation_levels
if (is.null(annotation_levels)) {
  predicted_cols <- grep("^predicted\\.", colnames(query[[]]), value = TRUE)
  annotation_levels <- sub(
    "^predicted\\.", "",
    predicted_cols[!grepl("\\.score$", predicted_cols)]
  )
}
level_names <- gsub(".", "_", annotation_levels, fixed = TRUE)
obs_predictions <- par$output_obs_predictions
if (is.null(obs_predictions)) {
  obs_predictions <- paste0("azimuth_pred_", level_names)
}
obs_probability <- par$output_obs_probability
if (is.null(obs_probability)) {
  obs_probability <- paste0("azimuth_probability_", level_names)
}

# Mapping of the RunAzimuth() column names to the output column names
obs_mapping <- c(
  setNames(obs_predictions, paste0("predicted.", annotation_levels)),
  setNames(
    obs_probability,
    paste0("predicted.", annotation_levels, ".score")
  ),
  setNames(par$output_obs_mapping_score, "mapping.score")
)
if (anyDuplicated(obs_mapping)) {
  stop(
    "The output column names should be unique, found: ",
    paste(obs_mapping, collapse = ", ")
  )
}
missing_cols <- setdiff(names(obs_mapping), colnames(query[[]]))
if (length(missing_cols) > 0) {
  stop(
    "RunAzimuth() did not return the expected columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

cat("Adding predicted labels and mapping scores to the cell metadata\n")
for (col in obs_mapping) {
  warn_overwrite(col, colnames(obj[[]]), "Metadata column")
}
result <- query[[]][cells, names(obs_mapping), drop = FALSE]
colnames(result) <- obs_mapping
obj <- AddMetaData(obj, metadata = result)

# Cell embeddings are not guaranteed to be in the same order as the cells of
# the object after IntegrateEmbeddings(), so realign them to the input object.
cat("Adding reference-projected reductions\n")
reduction_mapping <- c(
  "ref.umap" = par$output_obsm_umap,
  "integrated_dr" = par$output_obsm_integrated
)
for (reduction_name in names(reduction_mapping)) {
  output_name <- reduction_mapping[[reduction_name]]
  warn_overwrite(output_name, Reductions(obj), "Reduction")
  reduction <- query[[reduction_name]]
  reduction <- CreateDimReducObject(
    embeddings = Embeddings(reduction)[cells, , drop = FALSE],
    loadings = Loadings(reduction),
    key = paste0(gsub("[^[:alnum:]]", "", output_name), "_"),
    assay = par$assay,
    misc = Misc(reduction)
  )
  obj[[output_name]] <- reduction
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
