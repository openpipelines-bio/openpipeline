import scanpy as sc
import mudata as mu
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import re
import sys

## VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered_pseudobulk.h5mu",
    "output": "deseq2_results.csv",
    "input_layer": None,
    "modality": "rna",
    "obs_cell_group": None,
    "design_formula": "~ cell_type + disease + treatment",
    "contrast_column": "treatment",
    "contrast_values": [
        "stim",
        "ctrl",
    ],  # [treatment, baseline] or [group1, group2, group3, ...]
    "filter_genes_min_samples": None,
    "p_adj_threshold": 0.05,
    "log2fc_threshold": 0.0,
    "filter_gene_patterns": ["MIR\\d+", "AL\\d+", "LINC\\d+", "AC\\d+", "AP\\d+"],
    "var_gene_name": "feature_name",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])

from setup_logger import setup_logger

logger = setup_logger()


def is_normalized(layer):
    if sp.issparse(layer):
        row_sums = np.array(layer.sum(axis=1)).flatten()
    else:
        row_sums = layer.sum(axis=1)

    return np.allclose(row_sums, 1)


def parse_design_formula(design_formula):
    design_factors = re.findall(r"\b(?!~)\w+", design_formula)
    logger.info(
        f"Design formula: {design_formula}\nExtracted factors: {design_factors}"
    )
    return design_factors


def prepare_contrast_matrix(design_factors, contrast_column, metadata):
    # Validate required columns exist
    required_columns = (
        design_factors + [contrast_column]
        if contrast_column not in design_factors
        else design_factors
    )
    missing_columns = [col for col in required_columns if col not in metadata.columns]
    if missing_columns:
        raise ValueError(
            f"Missing required columns in metadata: {missing_columns}\n"
            f"Available metadata columns: {list(metadata.columns)}"
        )

    # Check contrast values exist
    contrast_values = par["contrast_values"]
    available_values = metadata[contrast_column].unique()
    missing_values = [val for val in contrast_values if val not in available_values]
    if missing_values:
        raise ValueError(
            f"Contrast values {missing_values} not found in {contrast_column}.\n"
            f"Available values: {available_values}"
        )

    if len(contrast_values) < 2:
        raise ValueError(f"Need at least 2 values for contrast, got: {contrast_values}")

    treatment, baseline = contrast_values
    contrast_tuple = (contrast_column, treatment, baseline)
    logger.info(
        f"Performing pairwise contrast: {contrast_column} {treatment} vs {baseline}"
    )
    return contrast_tuple


def prepare_counts_matrix(layer, mod):
    counts = pd.DataFrame(layer, columns=mod.var_names, index=mod.obs_names)
    counts = counts.astype(int)  # Ensure counts are integers (required for DESeq2)
    return counts


def filter_genes_by_pattern(counts, gene_pattern):
    pattern_string = "|".join(gene_pattern)
    before_filter = counts.shape[1]

    # Filter genes based on column names (gene names)
    genes_to_keep = list(
        ~pd.Series(counts.columns).str.contains(pattern_string, na=False)
    )
    counts = counts.loc[:, genes_to_keep]

    logger.info(
        f"Filtered out genes matching patterns {par['filter_gene_patterns']}:\n"
        f"Counts matrix shape before bene pattern filtering: {before_filter}\n"
        f"Counts matrix shape after gene pattern filtering: {counts.shape}"
    )

    return counts


def deseq2_analysis(counts, metadata, contrast_tuple, design_formula):
    # Creating DESeq2 dataset
    logger.info("Creating DESeq2 dataset")
    adata = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design_formula,
    )

    # Filtering genes based on presence across samples
    sample_count = (
        par["filter_genes_min_samples"] if par["filter_genes_min_samples"] else 1
    )
    logger.info(
        f"Filtering genes by counts: removing genes that are present in less than {sample_count} samples"
    )
    sc.pp.filter_genes(adata, min_cells=sample_count)

    # Run DESeq2 analysis
    logger.info("Running DESeq2 analysis")
    adata.deseq2()

    # Perform statistical test
    logger.info("Performing statistical test")
    stat_res = DeseqStats(adata, contrast=contrast_tuple)
    stat_res.summary()
    results = stat_res.results_df.reset_index()

    # Sort by log2FoldChange
    results = results.sort_values("log2FoldChange", ascending=False)

    # Add additional statistics
    results["abs_log2FoldChange"] = results["log2FoldChange"].abs()
    results["significant"] = (results["padj"] < par["p_adj_threshold"]) & (
        results["log2FoldChange"].abs() > par["log2fc_threshold"]
    )

    # Filter results based on significance
    significant_genes = results[
        (results["padj"] < par["p_adj_threshold"])
        & (results["log2FoldChange"].abs() > par["log2fc_threshold"])
    ]
    logger.info(
        f"Significant genes (padj < {par['p_adj_threshold']}, |log2FC| > {par['log2fc_threshold']}): {len(significant_genes)}"
    )
    return results


def main():
    # Load pseudobulk data
    logger.info(f"Loading pseudobulk data from {par['input']}")
    mdata = mu.read_h5mu(par["input"])
    mod = mdata.mod[par["modality"]]
    metadata = mod.obs.copy()
    layer = mod.layers[par["input_layer"]] if par["input_layer"] else mod.X
    if is_normalized(layer):
        raise ValueError("Input layer must contain raw counts.")

    # Parse design formula to extract factors
    logger.info("Preparing design formula")
    design_factors = parse_design_formula(par["design_formula"])

    # Preparing contrast matrix
    logger.info("Preparing contrast matrix")
    contrast_tuple = prepare_contrast_matrix(
        design_factors, par["contrast_column"], metadata
    )

    # Prepare counts matrix
    logger.info("Preparing counts matrix for DESeq2")
    counts = prepare_counts_matrix(layer, mod)

    # Filter out unwanted gene patterns if specified (before DESeq2 analysis)
    if par["filter_gene_patterns"]:
        logger.info("Filtering genes based on gene patterns")
        counts = filter_genes_by_pattern(counts, par["filter_gene_patterns"])

    # Run DESeq2 analysis
    try:
        # Check if running per cell group or multi-factor
        if par["obs_cell_group"]:
            logger.info("Running DESeq2 analysis per cell group")
            # Remove cell group from design formula for per-cell-type analysis
            design_no_celltype = (
                par["design_formula"]
                .replace(f"+ {par['obs_cell_group']}", "")
                .replace(f"{par['obs_cell_group']} +", "")
            )

            all_results = []
            for cell_group in metadata[par["obs_cell_group"]].unique():
                # Subset data for this cell group
                cell_mask = metadata[par["obs_cell_group"]] == cell_group
                counts_subset = counts.loc[cell_mask, :]
                metadata_subset = metadata.loc[cell_mask, :]

                # Run DESeq2 for this cell group
                logger.info(f"Running DESeq2 analysis for cell group {cell_group}...")
                cell_results = deseq2_analysis(
                    counts_subset, metadata_subset, contrast_tuple, design_no_celltype
                )

                # Add cell_group column to results
                cell_results["cell_group"] = cell_group
                all_results.append(cell_results)

            # Combine all results
            results = pd.concat(all_results, ignore_index=True)

        else:
            results = deseq2_analysis(
                counts, metadata, contrast_tuple, par["design_formula"]
            )

        # Log summary statistics
        upregulated = results[(results["log2FoldChange"] > 0) & results["significant"]]
        downregulated = results[
            (results["log2FoldChange"] < 0) & results["significant"]
        ]

        logger.info(
            "Summary:\n"
            f"  Total genes analyzed: {len(results)}\n"
            f"  Significant upregulated: {len(upregulated)}\n"
            f"  Significant downregulated: {len(downregulated)}\n"
        )

        # Save results
        logger.info(f"Saving results to {par['output']}")
        results.to_csv(par["output"], index=False)

    except Exception as e:
        logger.warning(
            "Check input data, design factors, and contrast specifications\n"
            f"Contrast column: {par['contrast_column']}\n"
            f"Contrast values: {par['contrast_values']}\n"
            f"Number of samples: {len(mod.obs)}\n"
            f"Number of genes: {len(mod.var)}\n"
        )

        raise e

    logger.info("DESeq2 analysis completed successfully")


if __name__ == "__main__":
    main()
