## VIASH START
par <- list(
  input = c(
    # "resources_test/bdrhap_5kjrt/processed/ABC_L1",
    # "resources_test/bdrhap_5kjrt/processed/ABC_L2",
    "resources_test/bdrhap_vdj/processed/VDJdemo"
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
met_groups <- map(par$input, function(dir) {
  file <- list.files(dir, pattern = "Metrics_Summary.csv$", full.names = TRUE)
  lines <- readr::read_lines(file)
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
        mutate(path = dir) %>%
        select(path, everything())
    }
  )
  names(group_dfs) <- group_titles

  group_dfs
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

cat("Reading in VDJ data\n")
# vdj_data <- lapply(par$input, function(dir) {
#   file <- list.files(dir, pattern = "_VDJ_perCell.csv$", full.names = TRUE)
#   readr::read_csv(
#     file,
#     comment = "#"
#   )
# })

cat("Constructing obs\n")
obs <- data.frame(
  row.names = unique_cells,
  sample_id = rep(par$id, length(unique_cells))
)

cat("Constructing var\n")
# var <- lapply(par$input, function(dir) {
#   file <- list.files(dir, pattern = "_Bioproduct_Stats.csv$", full.names = TRUE)
#   readr::read_csv(
#     file,
#     col_types = cols(.default = col_integer()),
#     comment = "#"
#   )
# })

var <- data.frame(
  row.names = unique_targets,
  feature_types = rep("Gene Expression", length(unique_targets))
  # todo: fill in feature types
)

cat("Constructing combined metrics summary\n")
process_qc <- function(name, aggr_fun) {
  raw <- map_df(met_groups, name)
  aggr <-
    if (nrow(raw) > 0) {
      tryCatch({
        aggr_fun(raw)
      }, error = function(e) {
        warning(paste0("Processing metrics table '", name, "' errored. Error message:\n", e$message))
        tibble()
      })
    } else {
      tibble()
    }
  comb <- bind_rows(raw, aggr)
  list(raw = raw, aggr = aggr, comb = comb)
}
# process sequencing quality
sequencing_quality <- process_qc("Sequencing Quality", function(raw) {
  raw %>%
    group_by(Library) %>%
    summarise(
      Total_Reads_in_FASTQ_ = sum(.data$Total_Reads_in_FASTQ),
      Pct_Read_Pair_Overlap = sum(.data$Total_Reads_in_FASTQ * .data$Pct_Read_Pair_Overlap) / Total_Reads_in_FASTQ_,
      Pct_Reads_Too_Short = sum(.data$Total_Reads_in_FASTQ * .data$Pct_Reads_Too_Short) / Total_Reads_in_FASTQ_,
      Pct_Reads_Low_Base_Quality = sum(.data$Total_Reads_in_FASTQ * .data$Pct_Reads_Low_Base_Quality) / Total_Reads_in_FASTQ_,
      Pct_Reads_High_SNF = sum(.data$Total_Reads_in_FASTQ * .data$Pct_Reads_High_SNF) / Total_Reads_in_FASTQ_,
      Pct_Reads_Filtered_Out = sum(.data$Total_Reads_in_FASTQ * .data$Pct_Reads_Filtered_Out) / Total_Reads_in_FASTQ_,
      Total_Reads_After_Quality_Filtering = sum(.data$Total_Reads_After_Quality_Filtering)
    ) %>%
    select(everything(), Library) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# process library quality
library_quality <- process_qc("Library Quality", function(raw) {
  raw %>%
    group_by(Library) %>%
    summarise(
      Total_Filtered_Reads_ = sum(.data$Total_Filtered_Reads),
      Pct_Contaminating_PhiX_Reads_in_Filtered_R2 = sum(.data$Total_Filtered_Reads * .data$Pct_Contaminating_PhiX_Reads_in_Filtered_R2) / Total_Filtered_Reads_,
      Pct_Q30_Bases_in_Filtered_R2 = sum(.data$Total_Filtered_Reads * .data$Pct_Q30_Bases_in_Filtered_R2) / Total_Filtered_Reads_,
      Pct_Assigned_to_Cell_Labels = sum(.data$Total_Filtered_Reads * .data$Pct_Assigned_to_Cell_Labels) / Total_Filtered_Reads_,
      Pct_Cellular_Reads_Aligned_Uniquely = sum(.data$Total_Filtered_Reads * .data$Pct_Cellular_Reads_Aligned_Uniquely) / Total_Filtered_Reads_
    ) %>%
    select(everything(), Library) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# process alignment categories
alignment_categories <- process_qc("Alignment Categories", function(raw) {
  raw %>%
    group_by(Library) %>%
    summarise(
      Cellular_Reads_ = sum(.data$Cellular_Reads),
      Annotated_Transcriptome_Pct = sum(.data$Cellular_Reads * .data$Annotated_Transcriptome_Pct) / Cellular_Reads_,
      Introns_Pct = sum(.data$Cellular_Reads * .data$Introns_Pct) / Cellular_Reads_,
      Intergenic_Regions_Pct = sum(.data$Cellular_Reads * .data$Intergenic_Regions_Pct) / Cellular_Reads_,
      Antisense_Pct = sum(.data$Cellular_Reads * .data$Antisense_Pct) / Cellular_Reads_,
      Not_Unique_Pct = sum(.data$Cellular_Reads * .data$Not_Unique_Pct) / Cellular_Reads_,
      Ambiguous_Pct = sum(.data$Cellular_Reads * .data$Ambiguous_Pct) / Cellular_Reads_,
      AbSeq_Pct = sum(.data$Cellular_Reads * .data$AbSeq_Pct) / Cellular_Reads_,
      Sample_Tag_Pct = sum(.data$Cellular_Reads * .data$Sample_Tag_Pct) / Cellular_Reads_,
      Unaligned_Pct = sum(.data$Cellular_Reads * .data$Unaligned_Pct) / Cellular_Reads_
    ) %>%
    select(everything(), Library) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# process reads and molecules
reads_and_molecules <- process_qc("Reads and Molecules", function(raw) {
  raw %>%
    group_by(Bioproduct_Type) %>%
    summarise(
      Aligned_Reads_By_Type_ = sum(.data$Aligned_Reads_By_Type),
      Total_Raw_Molecules_ = sum(.data$Total_Raw_Molecules),
      Total_RSEC_Molecules_ = sum(.data$Total_RSEC_Molecules),
      Total_DBEC_Molecules_ = sum(.data$Total_DBEC_Molecules),
      Mean_Raw_Sequencing_Depth = sum(.data$Total_Raw_Molecules * .data$Mean_Raw_Sequencing_Depth) / Total_Raw_Molecules_,
      Mean_RSEC_Sequencing_Depth = sum(.data$Total_RSEC_Molecules * .data$Mean_RSEC_Sequencing_Depth) / Total_RSEC_Molecules_,
      Mean_DBEC_Sequencing_Depth = sum(.data$Total_DBEC_Molecules * .data$Mean_DBEC_Sequencing_Depth) / Total_DBEC_Molecules_,
      Sequencing_Saturation = NA_real_,
      Pct_Cellular_Reads_with_Amplicons_Retained_by_DBEC = NA_real_
    ) %>%
    select(everything(), Bioproduct_Type) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# process cells rsec
cells_rsec <- process_qc("Cells RSEC", function(raw) {
  raw %>%
    group_by(Bioproduct_Type) %>%
    summarise(
      Putative_Cell_Count_ = sum(.data$Putative_Cell_Count),
      Pct_Reads_from_Putative_Cells = sum(.data$Putative_Cell_Count * .data$Mean_Reads_per_Cell * .data$Pct_Reads_from_Putative_Cells) / sum(.data$Putative_Cell_Count * .data$Mean_Reads_per_Cell),
      Mean_Reads_per_Cell = sum(.data$Putative_Cell_Count * .data$Mean_Reads_per_Cell) / Putative_Cell_Count_,
      Mean_Molecules_per_Cell = sum(.data$Putative_Cell_Count * .data$Mean_Molecules_per_Cell) / Putative_Cell_Count_,
      Median_Molecules_per_Cell = NA_real_,
      Mean_Bioproducts_per_Cell = sum(.data$Putative_Cell_Count * .data$Mean_Bioproducts_per_Cell) / Putative_Cell_Count_,
      Median_Bioproducts_per_Cell = NA_real_,
      Total_Bioproducts_Detected = length(unique_targets)
    ) %>%
    select(everything(), Bioproduct_Type) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# process cells dbec
cells_dbec <- process_qc("Cells DBEC", function(raw) {
  raw %>%
    group_by(Bioproduct_Type) %>%
    summarise(
      Putative_Cell_Count_ = sum(.data$Putative_Cell_Count),
      Pct_Reads_from_Putative_Cells = sum(.data$Putative_Cell_Count * .data$Mean_Reads_per_Cell * .data$Pct_Reads_from_Putative_Cells) / sum(.data$Putative_Cell_Count * .data$Mean_Reads_per_Cell),
      Mean_Reads_per_Cell = sum(.data$Putative_Cell_Count * .data$Mean_Reads_per_Cell) / Putative_Cell_Count_,
      Mean_Molecules_per_Cell = sum(.data$Putative_Cell_Count * .data$Mean_Molecules_per_Cell) / Putative_Cell_Count_,
      Median_Molecules_per_Cell = NA_real_,
      Mean_Bioproducts_per_Cell = sum(.data$Putative_Cell_Count * .data$Mean_Bioproducts_per_Cell) / Putative_Cell_Count_,
      Median_Bioproducts_per_Cell = NA_real_,
      Total_Bioproducts_Detected = length(unique_targets)
    ) %>%
    select(everything(), Bioproduct_Type) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# process error correction level
error_correction_level <- process_qc("Error Correction Level", function(raw) {
  raw %>%
    group_by(Bioproduct_Type) %>%
    summarise(
      Number_of_DBEC_and_RSEC_Corrected = NA_real_, # should this be a max or a sum?
      Number_of_RSEC_Corrected = NA_real_,
      Number_in_Panel = NA_real_
    ) %>%
    select(everything(), Bioproduct_Type) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# process vdj
vdj <- process_qc("VDJ", function(raw) {
  raw %>%
    group_by(Chain_Category) %>%
    summarise(
      Reads_Cellular_Aligned_to_VDJ_ = sum(.data$Reads_Cellular_Aligned_to_VDJ),
      Reads_Contig_Assembled_ = sum(.data$Reads_Contig_Assembled),
      Reads_VDJ_Annotated_ = sum(.data$Reads_VDJ_Annotated),
      Reads_Putative_ = sum(.data$Reads_Putative),
      Reads_Corrected_ = sum(.data$Reads_Corrected),
      Pct_Reads_Corrected = sum(.data$Pct_Reads_Corrected * .data$Reads_Contig_Assembled) / Reads_Contig_Assembled_,
      Mean_Reads_Corrected_per_Putative_Cell = mean(.data$Mean_Reads_Corrected_per_Putative_Cell), # assuming denominator is the same
      Molecules_VDJ_Annotated_ = sum(.data$Molecules_VDJ_Annotated),
      Molecules_Corrected_ = sum(.data$Molecules_Corrected),
      Mean_Molecules_Corrected_per_Putative_Cell = mean(.data$Mean_Molecules_Corrected_per_Putative_Cell), # assuming denominator is the same
      Dominant_Contigs_Mean_Nucleotide_Length = mean(.data$Dominant_Contigs_Mean_Nucleotide_Length), # assuming denominator is the same
      Dominant_Contigs_Pct_Full_Length = mean(.data$Dominant_Contigs_Pct_Full_Length), # assuming denominator is the same
      Dominant_Contigs_Pct_With_CDR3 = mean(.data$Dominant_Contigs_Pct_With_CDR3) # assuming denominator is the same
    ) %>%
    select(everything(), Chain_Category) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# process sample tags
sample_tags <- process_qc("Sample_Tags", function(raw) {
  raw %>%
    group_by(Chain_Category) %>%
    summarise(
      Sample_Tag_Filtered_Reads_ = sum(.data$Sample_Tag_Filtered_Reads),
      ST_Pct_Reads_from_Putative_Cells = NA_real_
    ) %>%
    select(everything(), Chain_Category) %>%
    rename_all(function(name) gsub("_$", "", name))
})

# combine mapping qcs
mapping_qc <- list(
  sequencing_quality = sequencing_quality,
  library_quality = library_quality,
  alignment_categories = alignment_categories,
  reads_and_molecules = reads_and_molecules,
  cells_rsec = cells_rsec,
  cells_dbec = cells_dbec,
  error_correction_level = error_correction_level,
  vdj = vdj,
  sample_tags = sample_tags
)
mapping_qc_per_file <- map(mapping_qc, "raw")
mapping_qc_combined <- map(mapping_qc, "comb")

cat("Constructing MuData object\n")
new_h5ad <- anndata::AnnData(
  X = new_counts,
  obs = obs,
  var = var,
  uns = list(
    mapping_qc_per_file = mapping_qc_per_file,
    mapping_qc_combined = mapping_qc_combined
  )
)

new_h5mu <- mudata$MuData(
  list(rna = new_h5ad)
)

cat("Writing to h5mu file\n")
new_h5mu$write(par$output)
