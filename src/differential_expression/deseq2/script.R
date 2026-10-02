library(DESeq2)
library(anndataR)
library(hdf5r)

rhdf5::h5disableFileLocking()

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
  p_adj_threshold = 0.05,
  log2fc_threshold = 0.0,
  var_gene_names = "feature_name",
  var_gene_symbol_column = NULL,
  export_normalized_counts = FALSE
)
meta <- list(resources_dir = "src/utils")
## VIASH END

cat("Starting DESeq2 analysis...\n")

h5mu_to_h5ad <- function(h5mu_path, modality_name) {
  tmp_path <- tempfile(fileext = ".h5ad")
  mod_location <- paste("mod", modality_name, sep = "/")
  h5src <- hdf5r::H5File$new(h5mu_path, "r")
  h5dest <- hdf5r::H5File$new(tmp_path, "w")
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

# Extract design factors from formula
parse_design_formula <- function(design_formula) {
  if (!grepl("^~\\s*\\w+(\\s*\\+\\s*\\w+)*$", design_formula)) {
    stop(
      sprintf(
        paste(
          "Invalid design formula: '%s'.",
          "Formula must start with '~' and contain factors separated by '+'.",
          "Example: '~ disease + treatment'",
          sep = "\n"
        ),
        design_formula
      )
    )
  }
  design_factors <- all.vars(as.formula(design_formula))
  cat("Design formula:", design_formula, "\n")
  cat("Extracted factors:", paste(design_factors, collapse = ", "), "\n")
  design_factors
}

# Validate and prepare contrast specifications
prepare_contrast_matrix <- function(
  design_factors, contrast_column, metadata
) {
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
    # Pairwise comparison: first value is the control group
    control_group <- contrast_values[1]
    comparison_group <- contrast_values[2]
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

# Create and configure DESeq2 dataset
create_deseq2_dataset <- function(
  counts, metadata, design_formula, gene_symbols = NULL
) {
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
  if (!is.null(gene_symbols)) {
    S4Vectors::mcols(dds)$gene_symbol <- gene_symbols
  }
  dds
}

# Extract DESeq2 results for each contrast from a fitted DESeqDataSet
deseq2_analysis <- function(dds, contrast_specs) {
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
    if (!is.null(S4Vectors::mcols(dds)$gene_symbol)) {
      results_df$gene_symbol <- S4Vectors::mcols(dds)$gene_symbol
    }
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

# Fit the DESeq2 model
fit_deseq2 <- function(dds) {
  cat("Running DESeq2 analysis\n")
  DESeq2::DESeq(dds)
}

# Write the per-sample table, normalized counts, variance-stabilized counts
# (VST, blind to the design) and run metadata of a fitted DESeqDataSet as
# "{file_prefix}_samples.csv", "{file_prefix}_normalized_counts.csv",
# "{file_prefix}_vst.csv" and "{file_prefix}_metadata.json".
export_normalized_counts <- function(
  dds, file_prefix, design_formula, contrast_specs, cell_group = NULL
) {
  cat("Exporting normalized and variance-stabilized counts\n")
  # A single contrast is a vector c(column, comparison, control)
  if (!is.list(contrast_specs)) {
    contrast_specs <- list(contrast_specs)
  }
  contrast_column <- contrast_specs[[1]][1]
  vst_function <- "vst"
  vst <- tryCatch(
    DESeq2::vst(dds, blind = TRUE),
    error = function(e) {
      # vst() fits the dispersion trend on a subset of 1000 genes, which fails
      # for small datasets; the full transformation does not have that limit
      cat(
        "vst() failed, using varianceStabilizingTransformation():",
        conditionMessage(e), "\n"
      )
      vst_function <<- "varianceStabilizingTransformation"
      DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
    }
  )

  # Samples: the metadata the design is built on, plus size factor and
  # library size (total raw counts)
  sample_columns <- intersect(
    unique(c(
      all.vars(as.formula(design_formula)),
      contrast_column,
      if (!is.null(cell_group)) par$obs_cell_group
    )),
    colnames(SummarizedExperiment::colData(dds))
  )
  samples <- data.frame(
    sample = colnames(dds),
    lapply(
      as.data.frame(SummarizedExperiment::colData(dds))[sample_columns],
      as.character
    ),
    size_factor = DESeq2::sizeFactors(dds),
    library_size = colSums(DESeq2::counts(dds)),
    check.names = FALSE
  )
  write.csv(samples, paste0(file_prefix, "_samples.csv"), row.names = FALSE)

  # Genes x samples tables
  gene_columns <- data.frame(gene_id = rownames(dds))
  if (!is.null(S4Vectors::mcols(dds)$gene_symbol)) {
    gene_columns$gene_symbol <- S4Vectors::mcols(dds)$gene_symbol
  }
  write_gene_table <- function(mat, suffix) {
    table <- data.frame(gene_columns, mat, check.names = FALSE)
    write.csv(
      table, paste0(file_prefix, "_", suffix, ".csv"),
      row.names = FALSE
    )
  }
  write_gene_table(DESeq2::counts(dds, normalized = TRUE), "normalized_counts")
  write_gene_table(SummarizedExperiment::assay(vst), "vst")

  # Run metadata
  metadata <- list(
    input = basename(par$input),
    modality = par$modality,
    input_layer = par$input_layer,
    cell_group = if (is.null(cell_group)) {
      NULL
    } else {
      list(
        column = par$obs_cell_group, value = as.character(cell_group)
      )
    },
    design_formula = design_formula,
    contrast_column = contrast_column,
    contrasts = lapply(contrast_specs, function(spec) {
      list(
        name = paste0(spec[2], "_vs_", spec[3]),
        comparison_group = spec[2],
        control_group = spec[3]
      )
    }),
    test = "Wald",
    p_adjust_method = "BH",
    p_adj_threshold = par$p_adj_threshold,
    log2fc_threshold = par$log2fc_threshold,
    lfc_shrinkage = FALSE,
    variance_stabilization = list(function_name = vst_function, blind = TRUE),
    n_samples = ncol(dds),
    n_genes = nrow(dds),
    var_gene_names = par$var_gene_names,
    var_gene_symbol_column = par$var_gene_symbol_column,
    versions = list(
      DESeq2 = as.character(utils::packageVersion("DESeq2")),
      R = paste(R.version$major, R.version$minor, sep = ".")
    )
  )
  jsonlite::write_json(
    metadata, paste0(file_prefix, "_metadata.json"),
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
}

# Save results and print summary statistics
save_results_and_log_summary <- function(
  results, output_file, cell_group = NULL
) {
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

  # Prepare analysis components
  cat("Preparing design formula\n")
  design_factors <- parse_design_formula(par$design_formula)

  cat("Preparing contrast matrix\n")
  contrast_specs <- prepare_contrast_matrix(
    design_factors, par$contrast_column, metadata
  )

  cat("Preparing counts matrix for DESeq2\n")
  var_names <- if (!is.null(par$var_gene_names)) {
    if (par$var_gene_names %in% colnames(mod$var)) {
      mod$var[[par$var_gene_names]]
    } else {
      stop(
        sprintf(
          "var_gene_names '%s' not found in mod$var columns: %s",
          par$var_gene_names,
          paste(colnames(mod$var), collapse = ", ")
        )
      )
    }
  } else {
    mod$var_names
  }
  gene_symbols <- NULL
  if (!is.null(par$var_gene_symbol_column)) {
    if (!par$var_gene_symbol_column %in% colnames(mod$var)) {
      stop(sprintf(
        "var_gene_symbol_column '%s' not found in mod$var columns: %s",
        par$var_gene_symbol_column, paste(colnames(mod$var), collapse = ", ")
      ))
    }
    gene_symbols <- as.character(mod$var[[par$var_gene_symbol_column]])
  }
  obs_names <- mod$obs_names
  counts <- prepare_counts_matrix(layer, var_names, obs_names)

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
        run_per_cell_group_analysis(
          counts, metadata, contrast_specs, gene_symbols
        )
      } else {
        run_overall_analysis(counts, metadata, contrast_specs, gene_symbols)
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
run_per_cell_group_analysis <- function(
  counts, metadata, contrast_specs, gene_symbols
) {
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
      counts_subset, metadata_subset, design_no_celltype, gene_symbols
    )
    dds <- fit_deseq2(dds)
    results <- deseq2_analysis(dds, contrast_specs)
    results[[par$obs_cell_group]] <- cell_group

    # Save results
    safe_name <- gsub("[/ \\(\\)]", "_", as.character(cell_group))
    safe_name <- gsub("_+", "_", safe_name)
    file_prefix <- file.path(
      par$output_dir,
      paste0(par$output_prefix, "_", safe_name)
    )
    save_results_and_log_summary(
      results, paste0(file_prefix, ".csv"), cell_group
    )
    if (par$export_normalized_counts) {
      export_normalized_counts(
        dds, file_prefix, design_no_celltype, contrast_specs, cell_group
      )
    }
  }
}

# Run overall analysis (all samples together)
run_overall_analysis <- function(
  counts, metadata, contrast_specs, gene_symbols
) {
  dds <- create_deseq2_dataset(
    counts, metadata, par$design_formula, gene_symbols
  )
  dds <- fit_deseq2(dds)
  results <- deseq2_analysis(dds, contrast_specs)

  file_prefix <- file.path(par$output_dir, par$output_prefix)
  save_results_and_log_summary(results, paste0(file_prefix, ".csv"))
  if (par$export_normalized_counts) {
    export_normalized_counts(
      dds, file_prefix, par$design_formula, contrast_specs
    )
  }
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}
