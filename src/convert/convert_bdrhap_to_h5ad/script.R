## VIASH START
par <- list(
  input = c(
    "output/64300535HPB1003_S10_rep1_lane2/",
    "output/64300535HPB1003_S10_rep1_lane1/"
  ),
  id = c(
    "64300535HPB1003_S10_rep1_lane2",
    "64300535HPB1003_S10_rep1_lane1"
  ),
  output = "/home/rcannood/test.h5ad"
)
## VIASH END

cat("Loading libraries\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)

read_metric <- function(file) {
  lines <- readr::read_lines(file)
  ix <- grep("^#[^#]", lines)
  header <- lines[ix+1] %>% strsplit(",") %>% unlist
  values <- lines[ix+2] %>% strsplit(",") %>% unlist
  ix2 <- !duplicated(header)
  new_csv <- paste0(
    header[ix2] %>% paste(collapse = ","), "\n",
    values[ix2] %>% paste(collapse = ","), "\n"
  )
  readr::read_csv(new_csv)
}

ids <- map_chr(seq_along(par$input), function(i) {
  if (!is.null(par$id) && length(par$input) == length(par$id)) {
    par$id[[i]]
  } else {
    paste0("sample", i)
  }
})

cat("Reading in metric summaries\n")
mets <- map_df(par$input, function(dir) {
  list.files(dir, pattern = "Metrics_Summary.csv$", full.names = TRUE) %>% read_metric()
})
mets$id <- ids

cat("Reading in count data\n")
counts <- lapply(par$input, function(dir) {
  list.files(dir, pattern = "_RSEC_MolsPerCell.csv$", full.names = TRUE) %>%
    readr::read_csv(
      col_types = cols(.default = col_integer()),
      comment = "#"
    ) %>%
    tibble::column_to_rownames("Cell_Index") %>%
    as.matrix %>%
    Matrix::Matrix(sparse = TRUE)
})

obs <- map_df(seq_along(counts), function(i) {
  cell_index = rownames(counts[[i]])
  data.frame(
    row.names = paste0("sample", i, "_", cell_index),
    id = rep(ids[[i]], length(cell_index))
  )
})

cat("Constructing var\n")
targets <- map(counts, colnames)
unique_targets <- sort(unique(unlist(targets)))
var <- data.frame(
  row.names = unique_targets,
  feature_types = rep("Gene Expression", length(unique_targets))
)

cat("Constructing counts\n")
new_counts <- map(seq_along(counts), function(i) {
  mat <- counts[[i]]
  if (is(mat, "RsparseMatrix")) {
    j <- mat@j+1
    jmap <- match(colnames(mat), unique_targets)
    newj <- jmap[j]
    Matrix::sparseMatrix(
      p = mat@p,
      j = newj,
      x = mat@x,
      repr = "R",
      dims = c(nrow(mat), length(unique_targets)),
      dimnames = list(rownames(mat), unique_targets)
    )
  } else if (is(mat, "CsparseMatrix")) {
    pmap <- cumsum(c(TRUE, unique_targets %in% colnames(mat)))
    newp <- mat@p[pmap]
    Matrix::sparseMatrix(
      p = newp,
      i = mat@i+1,
      x = mat@x,
      repr = "C",
      dims = c(nrow(mat), length(unique_targets)),
      dimnames = list(rownames(mat), unique_targets)
    )
  }
}) %>% do.call(rbind, .)

cat("Constructing metrics summary\n")
new_met <- tibble(
  Total_Reads_in_FASTQ = sum(mets$Total_Reads_in_FASTQ),
  Pct_Reads_Too_Short = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Reads_Too_Short) / Total_Reads_in_FASTQ,
  Pct_Reads_Low_Base_Quality = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Reads_Low_Base_Quality) / Total_Reads_in_FASTQ,
  Pct_Reads_High_SNF = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Reads_High_SNF) / Total_Reads_in_FASTQ,
  Pct_Reads_Filtered_Out = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Reads_Filtered_Out) / Total_Reads_in_FASTQ,
  Total_Reads_After_Quality_Filtering = sum(mets$Total_Reads_After_Quality_Filtering),
  Library = unique(mets$Library),
  Total_Filtered_Reads = sum(mets$Total_Filtered_Reads),
  Pct_Contaminating_PhiX_Reads_in_Filtered_R2 = sum(mets$Total_Filtered_Reads * mets$Pct_Contaminating_PhiX_Reads_in_Filtered_R2) / Total_Filtered_Reads,
  Pct_Q30_Bases_in_Filtered_R2 = sum(mets$Total_Filtered_Reads * mets$Pct_Q30_Bases_in_Filtered_R2) / Total_Filtered_Reads,
  Pct_Assigned_to_Cell_Labels = sum(mets$Total_Filtered_Reads * mets$Pct_Assigned_to_Cell_Labels) / Total_Filtered_Reads,
  Pct_Cellular_Reads_Aligned_Uniquely_to_Annotated_Transcriptome = sum(mets$Total_Filtered_Reads * mets$Pct_Cellular_Reads_Aligned_Uniquely_to_Annotated_Transcriptome) / Total_Filtered_Reads,
  Pct_Cellular_Reads_Aligned_Uniquely_to_Other_Genomic_Regions = sum(mets$Total_Filtered_Reads * mets$Pct_Cellular_Reads_Aligned_Uniquely_to_Other_Genomic_Regions) / Total_Filtered_Reads,
  Pct_Cellular_Reads_Aligned_Not_Unique = sum(mets$Total_Filtered_Reads * mets$Pct_Cellular_Reads_Aligned_Not_Unique) / Total_Filtered_Reads,
  Pct_Cellular_Reads_Unaligned = sum(mets$Total_Filtered_Reads * mets$Pct_Cellular_Reads_Unaligned) / Total_Filtered_Reads,
  Aligned_Reads_By_Type = sum(mets$Aligned_Reads_By_Type),
  Total_Raw_Molecules = sum(mets$Total_Raw_Molecules),
  Total_RSEC_Molecules = sum(mets$Total_RSEC_Molecules),
  Mean_Raw_Sequencing_Depth = sum(mets$Total_Raw_Molecules * mets$Mean_Raw_Sequencing_Depth) / Total_Raw_Molecules,
  Mean_RSEC_Sequencing_Depth = sum(mets$Total_RSEC_Molecules * mets$Mean_RSEC_Sequencing_Depth) / Total_RSEC_Molecules,
  Sequencing_Saturation = NA_real_,
  Target_Type = unique(mets$Target_Type),
  Putative_Cell_Count = sum(mets$Putative_Cell_Count),
  Pct_Reads_from_Putative_Cells = sum(mets$Putative_Cell_Count * mets$Mean_Reads_per_Cell * mets$Pct_Reads_from_Putative_Cells) / sum(mets$Putative_Cell_Count * mets$Mean_Reads_per_Cell),
  Mean_Reads_per_Cell = sum(mets$Putative_Cell_Count * mets$Mean_Reads_per_Cell) / Putative_Cell_Count,
  Mean_Molecules_per_Cell = sum(mets$Putative_Cell_Count * mets$Mean_Molecules_per_Cell) / Putative_Cell_Count,
  Median_Molecules_per_Cell = median(Matrix::rowSums(new_counts)),
  Mean_Targets_per_Cell = sum(mets$Putative_Cell_Count * mets$Mean_Targets_per_Cell) / Putative_Cell_Count,
  Median_Targets_per_Cell = median(Matrix::rowSums(new_counts > 0)),
  Total_Targets_Detected = length(unique_targets)
)

cat("Constructing anndata object\n")
new_h5ad <- anndata::AnnData(
  X = new_counts,
  obs = obs,
  var = var,
  uns = list(
    metrics_summary = new_met,
    metrics_per_file = mets %>% select(id, everything())
  )
)

cat("Storing raw\n")
new_h5ad$raw <- new_h5ad

cat("Writing to h5ad file\n")
new_h5ad$write_h5ad(par$output, compression = "gzip")
