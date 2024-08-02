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

cat("Reading in h5mu\n")
ori_h5mu_path <- list.files(par$input, pattern = "h5mu$", full.names = TRUE)
assert_that(
  length(ori_h5mu_path) == 1, 
  msg = paste0("Exactly one .h5mu should be found, found ", length(ori_h5mu_path), " h5mu files instead.")
)
ori_h5mu <- mudata$read_h5mu(ori_h5mu_path)

cat("Reading in count data")
counts_folder <- list.files(par$input, pattern = "_RSEC_MolsPerCell_MEX", full.names = TRUE, recursive = TRUE, include.dirs = TRUE)
assert_that(
length(counts_folder) == 1, 
msg = paste0("Exactly one *_RSEC_MolsPerCell_Unfiltered_MEX folder should be found, found ", length(counts_folder), " folders instead.")
)
count_files <- unzip(counts_folder, list = TRUE)


cat("Reading in feature file")
feature_file <- count_files$Name[str_detect(count_files$Name, "features.tsv.gz")]
assert_that(
    length(feature_file) == 1, 
    msg = paste0("Exactly one feature file should be found, found ", length(feature_file), " files instead.")
    )
tmp_file = tempfile()
unzip(counts_folder, files=feature_file, exdir=tmp_file)
features <-
  readr::read_tsv(
    file.path(tmp_file, feature_file),
    col_names = c("index", "feature_id", "feature_types")
    )

cat("Reading in barcodes files")
barcodes_file <- count_files$Name[str_detect(count_files$Name, "barcodes")]
assert_that(
    length(barcodes_file) == 1, 
    msg = paste0("Exactly one feature file should be found, found ", length(barcodes_file), " files instead.")
    )
tmp_file = tempfile()
unzip(counts_folder, files=barcodes_file, exdir=tmp_file)
barcodes <-
  readr::read_tsv(
    file.path(tmp_file, barcodes_file),
    col_names = c("cell_index")
    )

### PROCESSING RNA MODALITY
rna_adata <- ori_h5mu$mod[["rna"]]

cat("Constructing obs\n")
library_id <- metric_dfs[["sequencing_quality"]]$Library
if (length(library_id) > 1) {
  library_id <- paste(library_id[library_id != "Combined_stats"], collapse = " & ")
}

assert_that(
    identical(as.character(barcodes$cell_index), rownames(rna_adata$X)), 
    msg = paste0("The index of the rna modality X layer should be equal to the cell index column of the barcodes file.")
    )

rna_adata$obs$cell_id = rownames(rna_adata$X)
rna_adata$obs$run_id = par$id
rna_adata$obs$library_id = library_id
rna_adata$obs$sample_id = library_id

cat("Constructing var\n")

assert_that(
    all.equal(features$feature_id, rownames(rna_adata$var)), 
    msg = paste0("The index of the rna modality var index should be equal to the feature_id column of the features file.")
    )

rna_adata$var$gene_ids = par$id
rna_adata$var$gene_name = rownames(rna_adata$var)
rna_adata$var$feature_types = features$feature_types

cat("Constructing uns\n")
names(metric_dfs) <- paste0("mapping_qc_", names(metric_dfs))
rna_adata$uns <- c(metric_dfs)

#TODO: PROCESS OTHER MODALITIES

cat("Constructing MuData object\n")
modalities <-
  list(
    rna = rna_adata
  )

ori_h5mu <- mudata$MuData(modalities[!sapply(modalities, is.null)])

cat("Writing to h5mu file\n")
ori_h5mu$write(par$output, compression=par$output_compression)
