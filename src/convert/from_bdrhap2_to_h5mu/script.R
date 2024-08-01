cat("Loading libraries\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
library(assertthat)
library(Matrix)
library(unzip)
mudata <- reticulate::import("mudata")

## VIASH START
par <- list(
  id = "foo",
  input = "output_large",
  output = "bdrhap2_test.h5mu"
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

counts_folder <- list.files(par$input, pattern = "_RSEC_MolsPerCell_MEX", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
assert_that(
length(counts_folder) == 1, 
msg = paste0("Exactly one *_RSEC_MolsPerCell_Unfiltered_MEX folder should be found, found ", length(counts_folder), " folders instead.")
)

matrix_file <- count_files$Name[str_detect(count_files$Name, "matrix.mtx.gz")]
assert_that(
    length(matrix_file) == 1, 
    msg = paste0("Exactly one matrix file should be found, found ", length(matrix_file), " files instead.")
    )
tmp_file = tempfile()
unzip(counts_folder, files=matrix_file, exdir=tmp_file)
sparse_counts <- readMM(file.path(tmp_file, matrix_file))
unlink(tmp_file, recursive = TRUE)
counts <- t(
    as.data.frame(as.matrix(sparse_counts))
    )

feature_file <- count_files$Name[str_detect(count_files$Name, "features.tsv.gz")]
assert_that(
    length(matrix_file) == 1, 
    msg = paste0("Exactly one feature file should be found, found ", length(feature_file), " files instead.")
    )
tmp_file = tempfile()
unzip(counts_folder, files=feature_file, exdir=tmp_file)
features <-
  readr::read_tsv(
    file.path(tmp_file, feature_file),
    col_names = c("index", "feature_id", "feature_types")
    )

colnames(counts) <- as.character(features$feature_id)

cat("Constructing obs\n")
library_id <- metric_dfs[["sequencing_quality"]]$Library
if (length(library_id) > 1) {
  library_id <- paste(library_id[library_id != "Combined_stats"], collapse = " & ")
}

obs <- tibble(
  cell_id = rownames(counts),
  run_id = par$id,
  library_id = library_id,
  sample_id = library_id
)

cat("Constructing var\n")
# determine feature types of genes
var <- tibble(
    gene_ids = features$feature_id,
    gene_name = features$feature_id,
    feature_types = features$feature_types
    ) %>%
  as.data.frame() %>%
  column_to_rownames("gene_ids")

cat("Constructing uns\n")
names(metric_dfs) <- paste0("mapping_qc_", names(metric_dfs))
uns <- c(metric_dfs)

cat("Constructing RNA (&ABC?) AnnData")
adata <- anndata::AnnData(
  X = counts,
  obs = obs,
  var = var,
  uns = uns
)
adata_rna <- adata[, adata$var$feature_types != "Antibody Capture"]


cat("Constructing MuData object\n")
modalities <-
  list(
    rna = adata_rna
  )

mdata <- mudata$MuData(modalities[!sapply(modalities, is.null)])

cat("Writing to h5mu file\n")
mdata$write(par$output, compression=par$output_compression)

