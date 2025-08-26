library(DESeq2)
library(anndataR)
library(hdf5r)

## VIASH START
par <- list(
  input = paste0(
    "resources_test/annotation_test_data/",
    "TS_Blood_filtered_pseudobulk.h5mu"
  ),
  output_dir = "./test_deseq2_no_cellgroup/",
  output_prefix = "deseq2_analysis",
  input_layer = NULL,
  modality = "rna",
  obs_cell_group = "cell_type",
  design_formula = "~ treatment",
  contrast_column = "treatment",
  contrast_values = c("ctrl", "stim"),
  filter_genes_min_samples = NULL,
  p_adj_threshold = 0.05,
  log2fc_threshold = 0.0,
  filter_gene_patterns = c(
    "MIR\\d+", "AL\\d+", "LINC\\d+", "AC\\d+", "AP\\d+"
  ),
  var_gene_names = "feature_name"
)
meta <- list(resources_dir = "src/utils")
## VIASH END

cat("Starting DESeq2 analysis...\n")

h5mu_to_h5ad <- function(h5mu_path, modality_name) {
  tmp_path <- tempfile(fileext = ".h5ad")
  mod_location <- paste("mod", modality_name, sep = "/")
  h5src <- hdf5r::H5File$new(h5mu_path, "r")
  h5dest <- hdf5r::H5File$new(tmp_path, "w")
  # Copy over the child objects and the child attributes from root
  # Root cannot be copied directly because it always exists and
  # copying does not allow overwriting.
  children <- hdf5r::list.objects(h5src,
    path = mod_location,
    full.names = FALSE, recursive = FALSE
  )
  for (child in children) {
    h5dest$obj_copy_from(
      h5src, paste(mod_location, child, sep = "/"),
      paste0("/", child)
    )
  }
  # Also copy the root attributes
  root_attrs <- hdf5r::h5attr_names(x = h5src)
  for (attr in root_attrs) {
    h5a <- h5src$attr_open(attr_name = attr)
    robj <- h5a$read()
    h5dest$create_attr_by_name(
      attr_name = attr,
      obj_name = ".",
      robj = robj,
      space = h5a$get_space(),
      dtype = h5a$get_type()
    )
  }
  h5src$close()
  h5dest$close()

  tmp_path
}

# Check if expression data is normalized (row sums =~ 1)
is_normalized <- function(layer) {
  row_sums <- if (is(layer, "sparseMatrix") || is(layer, "dgCMatrix")) {
    Matrix::rowSums(layer)
  } else {
    rowSums(layer)
  }

  all(abs(row_sums - 1) < 1e-6, na.rm = TRUE)
}

# Extract design factors from formula
parse_design_formula <- function(design_formula) {
  design_factors <- all.vars(as.formula(design_formula))
  cat("Design formula:", design_formula, "\n")
  cat("Extracted factors:", paste(design_factors, collapse = ", "), "\n")
  design_factors
}

# Validate and prepare contrast specifications
prepare_contrast_matrix <- function(
    design_factors, contrast_column, metadata) {
  # Validate required columns exist
  required_columns <- unique(c(design_factors, contrast_column))
  missing_columns <- setdiff(required_columns, colnames(metadata))

  if (length(missing_columns) > 0) {
    stop(sprintf(
      paste(
        "Missing required columns in metadata: %s\n",
        "Available metadata columns: %s"
      ),
      paste(missing_columns, collapse = ", "),
      paste(colnames(metadata), collapse = ", ")
    ))
  }

  # Check contrast values exist
  contrast_values <- par$contrast_values
  available_values <- unique(metadata[[contrast_column]])
  missing_values <- setdiff(contrast_values, available_values)

  if (length(missing_values) > 0) {
    stop(sprintf(
      paste("Contrast values %s not found in %s.",
        "Available values: %s",
        sep = "\n"
      ),
      paste(missing_values, collapse = ", "),
      contrast_column,
      paste(available_values, collapse = ", ")
    ))
  }

  # Handle different contrast scenarios
  if (length(contrast_values) == 2) {
    # Pairwise comparison
    comparison_group <- contrast_values[1]
    control_group <- contrast_values[2]
    contrast_spec <- c(contrast_column, comparison_group, control_group)
    cat(
      "Performing pairwise contrast:", contrast_column,
      comparison_group, "vs", control_group, "\n"
    )
    contrast_spec
  } else if (length(contrast_values) > 2) {
    # Multiple comparisons against first value (control_group)
    control_group <- contrast_values[1]
    contrast_specs <- list()
    for (i in 2:length(contrast_values)) {
      comparison_group <- contrast_values[i]
      contrast_specs[[length(contrast_specs) + 1]] <-
        c(contrast_column, comparison_group, control_group)
    }
    cat(
      "Performing multiple contrasts against control_group '",
      control_group, "':",
      paste(sapply(contrast_specs, function(x) x[2]), collapse = ", "), "\n"
    )
    contrast_specs
  } else {
    stop(sprintf(
      "Need at least 2 values for contrast, got: %s",
      paste(contrast_values, collapse = ", ")
    ))
  }
}

# Convert expression matrix to counts data frame
prepare_counts_matrix <- function(layer, var_names, obs_names) {
  counts <- if (is(layer, "sparseMatrix") || is(layer, "dgCMatrix")) {
    as.matrix(layer)
  } else {
    layer
  }

  # Create properly named data frame (transpose for DESeq2 format)
  counts_df <- data.frame(counts)
  colnames(counts_df) <- var_names
  rownames(counts_df) <- obs_names

  # Ensure integer counts (required for DESeq2)
  counts_df[] <- lapply(counts_df, function(x) as.integer(round(x)))
  counts_df
}

# Filter genes based on regex patterns
filter_genes_by_pattern <- function(counts, gene_patterns) {
  if (is.null(gene_patterns) || length(gene_patterns) == 0) {
    return(counts)
  }

  pattern_string <- paste(gene_patterns, collapse = "|")
  before_filter <- ncol(counts)
  genes_to_keep <- !grepl(pattern_string, colnames(counts))
  counts_filtered <- counts[, genes_to_keep, drop = FALSE]

  cat(
    "Filtered out genes matching patterns:",
    paste(gene_patterns, collapse = ", "), "\n"
  )
  cat(
    "Matrix shape before filtering:", nrow(counts), "x", before_filter,
    "\n"
  )
  cat(
    "Matrix shape after filtering:", nrow(counts_filtered), "x",
    ncol(counts_filtered), "\n"
  )

  counts_filtered
}

# Create and configure DESeq2 dataset
create_deseq2_dataset <- function(counts, metadata, design_formula) {
  cat("Creating DESeq2 dataset\n")

  # Ensure matching samples between counts and metadata
  common_samples <- intersect(rownames(counts), rownames(metadata))
  if (length(common_samples) == 0) {
    stop("No common samples found between counts and metadata")
  }

  counts <- counts[common_samples, , drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  # Create DESeqDataSet (transpose for gene x sample format)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = t(counts),
    colData = metadata,
    design = as.formula(design_formula)
  )

  # Filtering genes based on presence across samples
  sample_count <- if (
    !is.null(par$filter_genes_min_samples)
  ) {
    par$filter_genes_min_samples
  } else {
    1
  }
  cat(
    "Filtering genes by counts: removing genes that are present in less than",
    sample_count, "samples\n"
  )

  keep <- rowSums(counts(dds) >= 1) >= sample_count
  dds <- dds[keep, ]

  dds
}

# Perform DESeq2 differential expression analysis
deseq2_analysis <- function(dds, contrast_specs) {
  cat("Running DESeq2 analysis\n")

  dds <- DESeq2::DESeq(dds)

  # Ensure contrast_specs is a list
  if (!is.list(contrast_specs)) {
    contrast_specs <- list(contrast_specs)
  }

  all_results <- lapply(seq_along(contrast_specs), function(i) {
    contrast_spec <- contrast_specs[[i]]
    cat(
      "Performing statistical test for contrast:",
      paste(contrast_spec, collapse = " "), "\n"
    )

    # Get DESeq2 results for this contrast
    res <- DESeq2::results(
      dds,
      contrast = contrast_spec, alpha = par$p_adj_threshold
    )

    # Convert to data frame and add metadata
    results_df <- as.data.frame(res)
    results_df$gene_id <- rownames(results_df)
    results_df$contrast <- paste0(contrast_spec[2], "_vs_", contrast_spec[3])
    results_df$comparison_group <- contrast_spec[2]
    results_df$control_group <- contrast_spec[3]
    results_df$abs_log2FoldChange <- abs(results_df$log2FoldChange)
    results_df$significant <- (
      results_df$padj < par$p_adj_threshold &
        !is.na(results_df$padj) &
        abs(results_df$log2FoldChange) > par$log2fc_threshold
    )

    # Sort by effect size
    results_df[order(results_df$log2FoldChange, decreasing = TRUE), ]
  })

  # Combine all results
  combined_results <- do.call(rbind, all_results)

  # Log summary statistics
  for (i in seq_along(contrast_specs)) {
    contrast_spec <- contrast_specs[[i]]
    contrast_name <- paste0(contrast_spec[2], "_vs_", contrast_spec[3])
    contrast_subset <- combined_results[
      combined_results$contrast == contrast_name,
    ]
    n_significant <- sum(contrast_subset$significant, na.rm = TRUE)
    cat("Contrast", contrast_name, ":", n_significant, "significant genes\n")
  }

  combined_results
}

# Save results and print summary statistics
save_results_and_log_summary <- function(
    results, output_file, cell_group = NULL) {
  group_text <- if (!is.null(cell_group)) paste(" for", cell_group) else ""
  cat("Saving results", group_text, "to", output_file, "\n")

  write.csv(results, output_file, row.names = FALSE)

  # Calculate summary statistics
  sig_results <- results[results$significant & !is.na(results$significant), ]
  upregulated <- sig_results[sig_results$log2FoldChange > 0, ]
  downregulated <- sig_results[sig_results$log2FoldChange < 0, ]

  cat("Summary", group_text, ":\n")
  cat("  Total genes analyzed:", nrow(results), "\n")
  cat("  Significant upregulated:", nrow(upregulated), "\n")
  cat("  Significant downregulated:", nrow(downregulated), "\n")
}

# Main analysis workflow
main <- function() {
  cat("Loading pseudobulk data from", par$input, "\n")

  # Load and prepare data
  h5ad_path <- h5mu_to_h5ad(par$input, par$modality)
  mod <- anndataR::read_h5ad(h5ad_path, as = "InMemoryAnnData")
  metadata <- as.data.frame(mod$obs)

  # Get expression matrix
  layer <- if (!is.null(par$input_layer)) {
    mod$layers[[par$input_layer]]
  } else {
    mod$X
  }

  if (is_normalized(layer)) {
    stop("Input layer must contain raw counts.")
  }

  # Prepare analysis components
  cat("Preparing design formula\n")
  design_factors <- parse_design_formula(par$design_formula)

  cat("Preparing contrast matrix\n")
  contrast_specs <- prepare_contrast_matrix(
    design_factors, par$contrast_column, metadata
  )

  cat("Preparing counts matrix for DESeq2\n")
  var_names <- if (
    !is.null(par$var_gene_names) && par$var_gene_names %in% colnames(mod$var)
  ) {
    mod$var[[par$var_gene_names]]
  } else {
    mod$var_names
  }
  obs_names <- mod$obs_names
  counts <- prepare_counts_matrix(layer, var_names, obs_names)

  # Apply gene filtering if specified
  if (
    !is.null(par$filter_gene_patterns) &&
      length(par$filter_gene_patterns) > 0
  ) {
    cat("Filtering genes based on patterns\n")
    counts <- filter_genes_by_pattern(counts, par$filter_gene_patterns)
  }

  # Ensure output directory exists
  if (!dir.exists(par$output_dir)) {
    dir.create(par$output_dir, recursive = TRUE)
  }

  # Run analysis (per cell group or overall)
  tryCatch(
    {
      if (
        !is.null(par$obs_cell_group) &&
          par$obs_cell_group %in% colnames(metadata)
      ) {
        run_per_cell_group_analysis(counts, metadata, contrast_specs)
      } else {
        run_overall_analysis(counts, metadata, contrast_specs)
      }
      cat("DESeq2 analysis completed successfully\n")
    },
    error = function(e) {
      cat("Error in analysis. Check input data and parameters:\n")
      cat("Contrast column:", par$contrast_column, "\n")
      cat("Contrast values:", paste(par$contrast_values, collapse = ", "), "\n")
      cat("Number of samples:", nrow(metadata), "\n")
      cat("Number of genes:", ncol(counts), "\n")
      stop(e)
    }
  )
}

# Run analysis per cell group
run_per_cell_group_analysis <- function(counts, metadata, contrast_specs) {
  cat("Running DESeq2 analysis per cell group\n")

  # Remove cell group from design formula
  design_no_celltype <- gsub(
    paste0("\\+\\s*", par$obs_cell_group), "", par$design_formula
  )
  design_no_celltype <- gsub(
    paste0(par$obs_cell_group, "\\s*\\+"), "", design_no_celltype
  )
  design_no_celltype <- gsub("\\s+", " ", design_no_celltype)

  cell_groups <- unique(metadata[[par$obs_cell_group]])

  for (cell_group in cell_groups) {
    cat("Processing cell group:", cell_group, "\n")

    # Subset data
    cell_mask <- metadata[[par$obs_cell_group]] == cell_group
    counts_subset <- counts[cell_mask, , drop = FALSE]
    metadata_subset <- metadata[cell_mask, , drop = FALSE]

    # Skip if insufficient samples
    if (nrow(counts_subset) < 2) {
      cat("Skipping cell group", cell_group, "- too few samples\n")
      next
    }

    # Run analysis
    dds <- create_deseq2_dataset(
      counts_subset, metadata_subset, design_no_celltype
    )
    results <- deseq2_analysis(dds, contrast_specs)
    results[[par$obs_cell_group]] <- cell_group

    # Save results
    safe_name <- gsub("[/ \\(\\)]", "_", as.character(cell_group))
    safe_name <- gsub("_+", "_", safe_name)
    output_file <- file.path(
      par$output_dir,
      paste0(par$output_prefix, "_", safe_name, ".csv")
    )
    save_results_and_log_summary(results, output_file, cell_group)
  }
}

# Run overall analysis (all samples together)
run_overall_analysis <- function(counts, metadata, contrast_specs) {
  dds <- create_deseq2_dataset(counts, metadata, par$design_formula)
  results <- deseq2_analysis(dds, contrast_specs)

  output_file <- file.path(
    par$output_dir,
    paste0(par$output_prefix, ".csv")
  )
  save_results_and_log_summary(results, output_file)
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
