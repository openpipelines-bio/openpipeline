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
    "obs_cell_group": "cell_type",
    "design_formula": "~ disease + treatment",
    "contrast_column": "treatment",
    "contrast_values": [
        "ctrl",
        "stim",
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

    # Handle multi-level contrasts
    if len(contrast_values) == 2:
        # Pairwise comparison (current behavior)
        treatment, baseline = contrast_values
        contrast_tuple = (contrast_column, treatment, baseline)
        logger.info(
            f"Performing pairwise contrast: {contrast_column} {treatment} vs {baseline}"
        )
        return [contrast_tuple]

    elif len(contrast_values) > 2:
        # Multiple comparisons - all vs first (baseline)
        baseline = contrast_values[0]
        contrast_tuples = []
        for treatment in contrast_values[1:]:
            contrast_tuples.append((contrast_column, treatment, baseline))
        logger.info(
            f"Performing multiple contrasts against baseline '{baseline}': {[t[1] for t in contrast_tuples]}"
        )
        return contrast_tuples

    else:
        raise ValueError(f"Need at least 2 values for contrast, got: {contrast_values}")


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


def create_deseq2_dataset(counts, metadata, design_formula):
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

    return adata


def deseq2_analysis(adata, contrast_tuples):
    logger.info("Running DESeq2 analysis")
    adata.deseq2()

    all_results = []

    # Handle multiple contrasts
    if not isinstance(contrast_tuples, list):
        contrast_tuples = [contrast_tuples]

    for contrast_tuple in contrast_tuples:
        logger.info(f"Performing statistical test for contrast: {contrast_tuple}")
        stat_res = DeseqStats(adata, contrast=contrast_tuple)
        stat_res.summary()
        results = stat_res.results_df.reset_index()

        # Add contrast information
        results["contrast"] = f"{contrast_tuple[1]}_vs_{contrast_tuple[2]}"
        results["treatment"] = contrast_tuple[1]
        results["baseline"] = contrast_tuple[2]

        # Sort by log2FoldChange
        results = results.sort_values("log2FoldChange", ascending=False)

        # Add additional statistics
        results["abs_log2FoldChange"] = results["log2FoldChange"].abs()
        results["significant"] = (results["padj"] < par["p_adj_threshold"]) & (
            results["log2FoldChange"].abs() > par["log2fc_threshold"]
        )

        all_results.append(results)

    # Combine all contrast results
    combined_results = pd.concat(all_results, ignore_index=True)

    # Log summary per contrast
    for contrast_tuple in contrast_tuples:
        contrast_name = f"{contrast_tuple[1]}_vs_{contrast_tuple[2]}"
        contrast_results = combined_results[
            combined_results["contrast"] == contrast_name
        ]
        significant_genes = contrast_results[contrast_results["significant"]]
        logger.info(
            f"Contrast {contrast_name}: {len(significant_genes)} significant genes"
        )

    return combined_results


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
    contrast_tuples = prepare_contrast_matrix(
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
                cell_dds = create_deseq2_dataset(
                    counts_subset, metadata_subset, design_no_celltype
                )
                cell_results = deseq2_analysis(cell_dds, contrast_tuples)

                # Add cell_group column to results
                cell_results[par["obs_cell_group"]] = cell_group
                all_results.append(cell_results)

            # Combine all results
            results = pd.concat(all_results, ignore_index=True)

        else:
            dds = create_deseq2_dataset(counts, metadata, par["design_formula"])
            results = deseq2_analysis(dds, contrast_tuples)

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
