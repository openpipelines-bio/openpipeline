## VIASH START
par <- list(
  input = c(
    "resources_test/bdrhap_5kjrt/processed/WTA"
  ),
  id = "foo",
  output = "test.h5mu"
)
## VIASH END

cat("Loading libraries\n")
options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("anndata", quietly = TRUE)
requireNamespace("reticulate", quietly = TRUE)
mudata <- reticulate::import("mudata")

cat("Reading in metric summaries\n")
mets <- map_df(par$input, function(dir) {
  file <- list.files(dir, pattern = "Metrics_Summary.csv$", full.names = TRUE)
  lines <- readr::read_lines(file)
  ix <- grep("^#[^#]", lines)
  header <- lines[ix + 1] %>% strsplit(",") %>% unlist
  values <- lines[ix + 2] %>% strsplit(",") %>% unlist
  ix2 <- !duplicated(header)
  new_csv <- paste0(
    header[ix2] %>% paste(collapse = ","), "\n",
    values[ix2] %>% paste(collapse = ","), "\n"
  )
  met_file <- readr::read_csv(new_csv)
  met_file %>%
    mutate(
      path = dir
    ) %>%
    select(path, everything())
})

cat("Reading in count data\n")
counts <- lapply(par$input, function(dir) {
  file <- list.files(dir, pattern = "_RSEC_MolsPerCell.csv$", full.names = TRUE)
  readr::read_csv(
    file,
    col_types = cols(.default = col_integer()),
    comment = "#"
  ) %>%
    tibble::column_to_rownames("Cell_Index") %>%
    as.matrix %>%
    Matrix::Matrix(sparse = TRUE)
})

cells <- map(counts, rownames)
unique_cells <- unique(unlist(cells))
targets <- map(counts, colnames)
unique_targets <- unique(unlist(targets))

cat("Constructing combined counts data\n")
new_counts <-
  map_df(counts, function(cou) {
    tup <- as(cou, "dgTMatrix")
    tibble(
      i = match(rownames(tup), unique_cells)[tup@i + 1],
      j = match(colnames(tup), unique_targets)[tup@j + 1],
      x = tup@x
    )
  }) %>%
  group_by(i, j) %>%
  summarise(x = sum(x), .groups = "drop") %>% {
    Matrix::sparseMatrix(
      i = .$i,
      j = .$j,
      x = .$x,
      dims = c(length(unique_cells), length(unique_targets)),
      dimnames = list(unique_cells, unique_targets)
    )
  }

cat("Constructing obs\n")
obs <- data.frame(
  row.names = unique_cells,
  sample_id = rep(par$id, length(unique_cells))
)

cat("Constructing var\n")
var <- data.frame(
  row.names = unique_targets,
  feature_types = rep("Gene Expression", length(unique_targets))
)

cat("Constructing combined metrics summary\n")
new_met <- tibble(
  Total_Reads_in_FASTQ = sum(mets$Total_Reads_in_FASTQ),
  Pct_Read_Pair_Overlap = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Read_Pair_Overlap) / Total_Reads_in_FASTQ,
  Pct_Reads_Too_Short = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Reads_Too_Short) / Total_Reads_in_FASTQ,
  Pct_Reads_Low_Base_Quality = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Reads_Low_Base_Quality) / Total_Reads_in_FASTQ,
  Pct_Reads_High_SNF = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Reads_High_SNF) / Total_Reads_in_FASTQ,
  Pct_Reads_Filtered_Out = sum(mets$Total_Reads_in_FASTQ * mets$Pct_Reads_Filtered_Out) / Total_Reads_in_FASTQ,
  Total_Reads_After_Quality_Filtering = sum(mets$Total_Reads_After_Quality_Filtering),
  Library = paste0(unique(mets$Library), collapse = ","),
  Total_Filtered_Reads = sum(mets$Total_Filtered_Reads),
  Pct_Contaminating_PhiX_Reads_in_Filtered_R2 = sum(mets$Total_Filtered_Reads * mets$Pct_Contaminating_PhiX_Reads_in_Filtered_R2) / Total_Filtered_Reads,
  Pct_Q30_Bases_in_Filtered_R2 = sum(mets$Total_Filtered_Reads * mets$Pct_Q30_Bases_in_Filtered_R2) / Total_Filtered_Reads,
  Pct_Assigned_to_Cell_Labels = sum(mets$Total_Filtered_Reads * mets$Pct_Assigned_to_Cell_Labels) / Total_Filtered_Reads,
  Pct_Cellular_Reads_Aligned_Uniquely = sum(mets$Total_Filtered_Reads * mets$Pct_Cellular_Reads_Aligned_Uniquely) / Total_Filtered_Reads,
  Cellular_Reads = sum(mets$Cellular_Reads),
  Annotated_Transcriptome_Pct = sum(mets$Cellular_Reads * mets$Annotated_Transcriptome_Pct) / Cellular_Reads,
  Introns_Pct = sum(mets$Cellular_Reads * mets$Introns_Pct) / Cellular_Reads,
  Intergenic_Regions_Pct = sum(mets$Cellular_Reads * mets$Intergenic_Regions_Pct) / Cellular_Reads,
  Antisense_Pct = sum(mets$Cellular_Reads * mets$Antisense_Pct) / Cellular_Reads,
  Not_Unique_Pct = sum(mets$Cellular_Reads * mets$Not_Unique_Pct) / Cellular_Reads,
  Ambiguous_Pct = sum(mets$Cellular_Reads * mets$Ambiguous_Pct) / Cellular_Reads,
  AbSeq_Pct = sum(mets$Cellular_Reads * mets$AbSeq_Pct) / Cellular_Reads,
  Sample_Tag_Pct = sum(mets$Cellular_Reads * mets$Sample_Tag_Pct) / Cellular_Reads,
  Unaligned_Pct = sum(mets$Cellular_Reads * mets$Unaligned_Pct) / Cellular_Reads,
  Aligned_Reads_By_Type = sum(mets$Aligned_Reads_By_Type),
  Total_Raw_Molecules = sum(mets$Total_Raw_Molecules),
  Total_RSEC_Molecules = sum(mets$Total_RSEC_Molecules),
  Mean_Raw_Sequencing_Depth = sum(mets$Total_Raw_Molecules * mets$Mean_Raw_Sequencing_Depth) / Total_Raw_Molecules,
  Mean_RSEC_Sequencing_Depth = sum(mets$Total_RSEC_Molecules * mets$Mean_RSEC_Sequencing_Depth) / Total_RSEC_Molecules,
  Sequencing_Saturation = NA_real_,
  Bioproduct_Type = paste0(unique(mets$Bioproduct_Type), collapse = ","),
  Putative_Cell_Count = sum(mets$Putative_Cell_Count),
  Pct_Reads_from_Putative_Cells = sum(mets$Putative_Cell_Count * mets$Mean_Reads_per_Cell * mets$Pct_Reads_from_Putative_Cells) / sum(mets$Putative_Cell_Count * mets$Mean_Reads_per_Cell),
  Mean_Reads_per_Cell = sum(mets$Putative_Cell_Count * mets$Mean_Reads_per_Cell) / Putative_Cell_Count,
  Mean_Molecules_per_Cell = sum(mets$Putative_Cell_Count * mets$Mean_Molecules_per_Cell) / Putative_Cell_Count,
  Median_Molecules_per_Cell = median(Matrix::rowSums(new_counts)),
  Mean_Bioproducts_per_Cell = sum(mets$Putative_Cell_Count * mets$Mean_Bioproducts_per_Cell) / Putative_Cell_Count,
  Median_Bioproducts_per_Cell = median(Matrix::rowSums(new_counts > 0)),
  Total_Bioproducts_Detected = length(unique_targets)
) %>%
  as.data.frame()
  
cat("Constructing MuData object\n")
new_h5ad <- anndata::AnnData(
  X = new_counts,
  obs = obs,
  var = var,
  uns = list(
    sample_metrics = new_met
  )
)

new_h5mu <- mudata$MuData(
  list(rna = new_h5ad)
)

cat("Writing to h5mu file\n")
new_h5mu$write(par$output)
