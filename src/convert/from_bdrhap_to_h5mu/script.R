cat("Loading libraries\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
library(assertthat)
mudata <- reticulate::import("mudata")

## VIASH START
par <- list(
  id = "foo",
  input = "resources_test/bdrhap_5kjrt/processed/WTA.bd_rhapsody.output_raw",
  output = "test.h5mu"
)
## VIASH END

read_metrics <- function(file) {
  metric_lines <- readr::read_lines(file)
  metric_lines_no_header <- metric_lines[!grepl("^##", metric_lines)]

  # parse sub data frames
  group_title_regex <- "^#([^#]*)#"
  group_title_ix <- grep(group_title_regex, metric_lines_no_header)
  group_titles <- gsub(group_title_regex, "\\1", metric_lines_no_header[group_title_ix])
  group_ix_from <- group_title_ix+1
  group_ix_to <- c(group_title_ix[-1]-1, length(metric_lines_no_header))
  metric_dfs <- pmap(
    list(
      from = group_ix_from,
      to = group_ix_to
    ),
    function(from, to) {
      lines <- metric_lines_no_header[from:to]
      lines <- lines[lines != ""]
      readr::read_csv(paste0(lines, collapse = "\n")) %>%
        mutate(run_id = par$id) %>%
        select(run_id, everything())
    }
  )
  names(metric_dfs) <- gsub(" ", "_", tolower(group_titles))
  metric_dfs
}
cat("Reading in metric summaries\n")
metrics_file <- list.files(par$input, pattern = "_Metrics_Summary.csv$", full.names = TRUE)
assert_that(
  length(metrics_file) == 1, 
  msg = paste0("Exactly one *_Metrics_Summary.csv should be found, found ", length(metrics_file), " files instead.")
)
metric_dfs <- read_metrics(metrics_file)

cat("Reading in count data\n")
counts_file <- list.files(par$input, pattern = "_DBEC_MolsPerCell.csv$", full.names = TRUE)
if (length(counts_file) == 0) {
  cat("Warning: could not find DBEC file, looking for RSEC file instead.\n")
  counts_file <- list.files(par$input, pattern = "_RSEC_MolsPerCell.csv$", full.names = TRUE)
}
assert_that(
  length(counts_file) == 1,
  msg = paste0("Exactly one *_(RSEC|DBEC)_MolsPerCell.csv should be found, found ", length(counts_file), " files instead.")
)
counts <-
  readr::read_csv(
    counts_file,
    col_types = cols(.default = col_integer()),
    comment = "#"
  ) %>%
    tibble::column_to_rownames("Cell_Index") %>%
    as.matrix %>%
    Matrix::Matrix(sparse = TRUE)

# processing VDJ data
vdj_file <- list.files(par$input, pattern = "_VDJ_perCell.csv$", full.names = TRUE)
vdj_data <-
  if (length(vdj_file) == 1) {
    cat("Reading in VDJ data\n")
    readr::read_csv(
      vdj_file,
      comment = "#"
    )
  } else {
    NULL
  }

cat("Reading in VDJ metric summaries\n")
vdj_metrics_file <- list.files(par$input, pattern = "_VDJ_metrics.csv$", full.names = TRUE)
vdj_metric_dfs <-
  if (length(vdj_metrics_file) == 1) {
    read_metrics(vdj_metrics_file)
  } else {
    list()
  }

# processing SMK data
smk_file <- list.files(par$input, pattern = "_Sample_Tag_Calls.csv$", full.names = TRUE)
smk_calls <-
  if (length(smk_file) == 1) {
    cat("Processing sample tags\n")
    readr::read_csv(
      smk_file,
      comment = "#"
    )
  } else {
    NULL
  }
smk_metrics_file <- list.files(par$input, pattern = "_Sample_Tag_Metrics.csv$", full.names = TRUE)
smk_metrics <-
  if (length(smk_metrics_file) == 1) {
    readr::read_csv(
      smk_metrics_file,
      comment = "#"
    )
  } else {
    NULL
  }

cat("Constructing obs\n")
library_id <- metric_dfs[["sequencing_quality"]]$Library
if (length(library_id) > 1) {
  library_id <- paste(library_id[library_id != "Combined_stats"], collapse = " & ")
}

obs <- tibble(
  cell_id = rownames(counts),
  run_id = par$id,
  library_id = library_id
)

if (!is.null(smk_calls)) {
  obs <- left_join(
    obs,
    smk_calls %>% transmute(
      cell_id = as.character(Cell_Index),
      sample_tag = Sample_Tag,
      sample_id = Sample_Name
    ),
    by = "cell_id"
  )
} else {
  obs <- obs %>% mutate(sample_id = library_id)
}

obs <- obs %>%
  mutate(sample_id = ifelse(!is.na(sample_id), sample_id, run_id)) %>%
  as.data.frame() %>%
  column_to_rownames("cell_id")

cat("Constructing var\n")
# determine feature types of genes
var0 <- tryCatch({
  feature_types_file <- list.files(par$input, pattern = "feature_types.tsv$", full.names = TRUE)

  # abseq fasta reference has trailing info which apparently gets stripped off by the bd rhapsody pipeline
  readr::read_tsv(feature_types_file) %>%
    mutate(
      trimmed_feature_id = gsub(" .*", "", feature_id),
      i = match(feature_id, colnames(counts)),
      j = match(trimmed_feature_id, colnames(counts)),
      ij = ifelse(is.na(i), j, i),
      final_feature_id = ifelse(!is.na(i), feature_id, trimmed_feature_id)
    ) %>%
    filter(!is.na(ij)) %>%
    select(feature_id = final_feature_id, feature_type, reference_file)
}, error = function(e) {
  cat("Feature matching error: ", e$message, "\n", sep = "")
  tibble(
    feature_id = character()
  )
})

# in case the feature types are missing
missing_features <- tibble(
  feature_id = setdiff(colnames(counts), var0$feature_id),
  feature_type = "Gene Expression",
  reference_file = NA_character_,
  note = "Feature annotation file missing, assuming type is Gene Expression"
)
var1 <-
  if (nrow(missing_features) > 0) {
    cat("Feature annotation file missing, assuming type is Gene Expression\n")
    bind_rows(var0, missing_features) %>%
      slice(match(colnames(counts), feature_id))
  } else {
    var0
  }

# create var
var <- var1 %>%
  transmute(gene_ids = feature_id, gene_name = feature_id, feature_types = feature_type, reference_file) %>%
  as.data.frame() %>%
  column_to_rownames("gene_ids")

cat("Constructing uns\n")
names(metric_dfs) <- paste0("mapping_qc_", names(metric_dfs))
smk_metric_dfs <-
  if (!is.null(smk_metrics)) {
    list(mapping_qc_smk_metrics = smk_metrics)
  } else {
    NULL
  }
uns <- c(metric_dfs, smk_metric_dfs)

cat("Constructing RNA (&ABC?) AnnData")
adata <- anndata::AnnData(
  X = counts,
  obs = obs,
  var = var,
  uns = uns
)

adata_prot <- adata[, adata$var$feature_types == "Antibody Capture"]
if (ncol(adata_prot) == 0) {
  adata_prot <- NULL
}
adata_rna <- adata[, adata$var$feature_types != "Antibody Capture"]

adata_vdj <-
  if (!is.null(vdj_data)) {
    cat("Constructing VDJ AnnData\n")
    names(vdj_metric_dfs) <- paste0("mapping_qc_", names(vdj_metric_dfs))
    anndata::AnnData(
      obs = vdj_data,
      uns = vdj_metric_dfs,
      shape = c(nrow(vdj_data), 0L)
    )
  } else {
    NULL
  }

cat("Constructing MuData object\n")
modalities <-
  list(
    rna = adata_rna,
    prot = adata_prot,
    vdj = adata_vdj
  )
mdata <- mudata$MuData(modalities[!sapply(modalities, is.null)])

cat("Writing to h5mu file\n")
mdata$write(par$output, compression="gzip")
