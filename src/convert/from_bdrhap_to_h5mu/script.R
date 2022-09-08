cat("Loading libraries\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
mudata <- reticulate::import("mudata")

## VIASH START
par <- list(
  id = "foo",
  # input = "resources_test/bdrhap_5kjrt/processed/ABC",
  input = "resources_test/bdrhap_vdj/processed/VDJdemo",
  output = "test.h5mu"
)
## VIASH END

cat("Reading in metric summaries\n")
metrics_file <- list.files(par$input, pattern = "Metrics_Summary.csv$", full.names = TRUE)
lines <- readr::read_lines(metrics_file)
lines_no_header <- lines[!grepl("^##", lines)]

# parse sub data frames
group_title_regex <- "^#([^#]*)#"
group_title_ix <- grep(group_title_regex, lines_no_header)
group_titles <- gsub(group_title_regex, "\\1", lines_no_header[group_title_ix])
group_ix_from <- group_title_ix+1
group_ix_to <- c(group_title_ix[-1]-1, length(lines_no_header))
group_dfs <- pmap(
  list(
    from = group_ix_from,
    to = group_ix_to
  ),
  function(from, to) {
    lines <- lines_no_header[from:to]
    lines <- lines[lines != ""]
    readr::read_csv(paste0(lines, collapse = "\n")) %>%
      mutate(sample_id = par$id) %>%
      select(sample_id, everything())
  }
)
names(group_dfs) <- group_titles


cat("Reading in count data\n")
counts_file <- list.files(par$input, pattern = "_RSEC_MolsPerCell.csv$", full.names = TRUE)
counts <-
  readr::read_csv(
    counts_file,
    col_types = cols(.default = col_integer()),
    comment = "#"
  ) %>%
    tibble::column_to_rownames("Cell_Index") %>%
    as.matrix %>%
    Matrix::Matrix(sparse = TRUE)

cat("Reading in VDJ data\n")
vdj_file <- list.files(par$input, pattern = "_VDJ_perCell.csv$", full.names = TRUE)
vdj_counts <-
  if (length(vdj_file) == 1) {
    readr::read_csv(
      file,
      comment = "#"
    )
  } else {
    NULL
  }

cat("Constructing obs\n")
obs <- data.frame(
  row.names = rownames(counts),
  sample_id = rep(par$id, nrow(counts))
)

# todo: feature type of ABC should not be GEX!
cat("Constructing var\n")
bioproduct_file <- list.files(par$input, pattern = "_Bioproduct_Stats.csv$", full.names = TRUE)
rna_var <- readr::read_csv(
  bioproduct_file,
  comment = "#"
) %>%
  mutate(feature_types = rep("Gene Expression", nrow(counts)))

cat("Constructing MuData object\n")
rna_h5ad <- anndata::AnnData(
  X = counts,
  obs = obs,
  var = var,
  uns = list(
    mapping_qc_sequencing_quality = group_dfs[["Sequencing Quality"]],
    mapping_qc_library_quality = group_dfs[["Library Quality"]],
    mapping_qc_alignment_categories = group_dfs[["Alignment Categories"]],
    mapping_qc_reads_and_molecules = group_dfs[["Reads and Molecules"]],
    mapping_qc_cells_rsec = group_dfs[["Cells RSEC"]],
    mapping_qc_cells_dbec = group_dfs[["Cells DBEC"]],
    mapping_qc_error_correction_level = group_dfs[["Error Correction Level"]],
    mapping_qc_vdj = group_dfs[["VDJ"]],
    mapping_qc_sample_tags = group_dfs[["Sample_Tags"]]
  )
)

new_h5mu <- mudata$MuData(
  list(rna = rna_h5ad)
)

cat("Writing to h5mu file\n")
new_h5mu$write(par$output)
