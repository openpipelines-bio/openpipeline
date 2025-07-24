import scanpy as sc
import mudata as mu
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import re

## VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered_annotated_pseudobulk.h5mu",
    "output": "deseq2_results.csv",
    "input_layer": None,
    "design_formula": "~ cell_type + disease + treatment",
    "contrast_column": "treatment",
    "contrast_values": ["stim", "ctrl"],  # [treatment, baseline] or [group1, group2, group3, ...]
    "filter_genes_min_samples": None,
    "padj_threshold": 0.05,
    "log2fc_threshold": 0.0,
    "filter_gene_patterns": ["MIR\\d+", "AL\\d+", "LINC\\d+", "AC\\d+", "AP\\d+"],
    "var_gene_name" : "feature_name",
    "per_cell_type": False,
    "cell_type_column": "cell_type",
}
meta = {"resources_dir": "src/utils"}

## VIASH END

def is_normalized(layer):
    if sp.issparse(layer):
        row_sums = np.array(layer.sum(axis=1)).flatten()
    else:
        row_sums = layer.sum(axis=1)

    return np.allclose(row_sums, 1)


import sys
sys.path.append(meta["resources_dir"])

from setup_logger import setup_logger

logger = setup_logger()


# Load pseudobulk data
logger.info(f"Loading pseudobulk data from {par['input']}")
mdata = mu.read_h5mu(par["input"])
mod = mdata.mod["rna"]

logger.info(f"Input pseudobulk data shape: {mod.shape}")
logger.info(f"Available metadata columns: {list(mod.obs.columns)}") 

# Parse design formula to extract factors
design_factors = re.findall(r'\b(?!~)\w+', par["design_formula"])
logger.info(f"Design formula: {par['design_formula']}")
logger.info(f"Extracted factors: {design_factors}")

# Validate required columns exist
contrast_column = par["contrast_column"]
required_columns = design_factors + [contrast_column] if contrast_column not in design_factors else design_factors
missing_columns = [col for col in required_columns if col not in mod.obs.columns]
if missing_columns:
    raise ValueError(f"Missing required columns in metadata: {missing_columns}")

# Prepare counts matrix
logger.info("Preparing counts matrix for DESeq2")

# Use specified layer if provided
if par["input_layer"]:
    mod.X = mod.layers[par["input_layer"]]

if is_normalized(mod.X):
    raise ValueError("Input layer must contain raw counts.")

counts = pd.DataFrame(
    mod.X, 
    columns=mod.var_names,
    index=mod.obs_names
)

# Ensure counts are integers (required for DESeq2)
counts = counts.astype(int)

metadata = mod.obs.copy()

# Filter out unwanted gene patterns if specified (before DESeq2 analysis)
if par["filter_gene_patterns"]:
    pattern_string = "|".join(par["filter_gene_patterns"])
    before_filter = counts.shape[1]
    
    # Filter genes based on column names (gene names)
    genes_to_keep = ~pd.Series(counts.columns).str.contains(pattern_string, na=False)
    counts_filtered = counts.loc[:, genes_to_keep]
    
    logger.info(f"Filtered out genes matching patterns {par['filter_gene_patterns']}: {before_filter} -> {counts_filtered.shape[1]}")
    counts = counts_filtered

logger.info(f"Counts matrix shape after gene pattern filtering: {counts.shape}")
logger.info(f"Design formula: {par['design_formula']}")

# Check contrast values exist
contrast_column = par["contrast_column"]
contrast_values = par["contrast_values"]
available_values = metadata[contrast_column].unique()
logger.info(f"Available values in {contrast_column}: {available_values}")

# Validate all contrast values exist
missing_values = [val for val in contrast_values if val not in available_values]
if missing_values:
    raise ValueError(f"Contrast values {missing_values} not found in {contrast_column}")

if len(contrast_values) < 2:
    raise ValueError(f"Need at least 2 values for contrast, got: {contrast_values}")

treatment, baseline = contrast_values
contrast_tuple = (contrast_column, treatment, baseline)
logger.info(f"Performing pairwise contrast: {contrast_column} {treatment} vs {baseline}")

# Create DESeq2 dataset
logger.info("Creating DESeq2 dataset")
try:
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=par["design_formula"],
    )
    
    # Optional gene filtering
    if par["filter_genes_min_samples"] is not None:
        logger.info(f"Filtering genes expressed in at least {par['filter_genes_min_samples']} samples")
        sc.pp.filter_genes(dds, min_cells=par["filter_genes_min_samples"])
    else:
        # Default: only filter genes with zero counts across all samples
        logger.info("Default filtering: removing genes with zero counts across all samples")
        sc.pp.filter_genes(dds, min_cells=1)
    
    logger.info(f"After filtering: {dds.shape}")
    
    # Check if running per cell type or multi-factor
    if par["per_cell_type"]:
        logger.info("Running DESeq2 analysis per cell type")
        
        # Remove cell_type from design formula for per-cell-type analysis
        design_no_celltype = par["design_formula"].replace(f"+ {par['cell_type_column']}", "").replace(f"{par['cell_type_column']} +", "")
        
        all_results = []
        for cell_type in metadata[par["cell_type_column"]].unique():
            logger.info(f"Analyzing cell type: {cell_type}")
            
            # Subset data for this cell type
            cell_mask = metadata[par["cell_type_column"]] == cell_type
            counts_subset = counts.loc[cell_mask, :]
            metadata_subset = metadata.loc[cell_mask, :]
            
            # Run DESeq2 for this cell type
            dds_cell = DeseqDataSet(
                counts=counts_subset,
                metadata=metadata_subset,
                design=design_no_celltype,
            )
            
            # Run DESeq2 analysis
            logger.info("Running DESeq2 analysis...")
            dds_cell.deseq2()
            
            # Perform statistical test
            stat_res_cell = DeseqStats(
                dds_cell, 
                contrast=contrast_tuple
            )
            stat_res_cell.summary()
            
            # Get results
            cell_results = stat_res_cell.results_df.reset_index()
            logger.info(f"DESeq2 analysis for cell type {cell_type} completed. Found {len(cell_results)} genes.")
            
            # Filter results based on significance
            significant_genes_cell = cell_results[
                (cell_results['padj'] < par["padj_threshold"]) & 
                (cell_results['log2FoldChange'].abs() > par["log2fc_threshold"])
            ]
            logger.info(f"Significant genes for cell type {cell_type} (padj < {par['padj_threshold']}, |log2FC| > {par['log2fc_threshold']}): {len(significant_genes_cell)}")
            
            # Sort by log2FoldChange
            cell_results = cell_results.sort_values('log2FoldChange', ascending=False)
            
            # Add additional statistics
            cell_results['abs_log2FoldChange'] = cell_results['log2FoldChange'].abs()
            cell_results['significant'] = (
                (cell_results['padj'] < par["padj_threshold"]) & 
                (cell_results['log2FoldChange'].abs() > par["log2fc_threshold"])
            )
            
            # Add cell_type column to results
            cell_results['cell_type'] = cell_type
            all_results.append(cell_results)
        
        # Combine all results
        results = pd.concat(all_results, ignore_index=True)
        
    else:
        logger.info("Running multi-factor DESeq2 analysis")
        logger.info("Running DESeq2 analysis...")
        dds.deseq2()
        
        # Perform statistical test
        stat_res = DeseqStats(
            dds, 
            contrast=contrast_tuple
        )
        stat_res.summary()
        
        # Get results
        results = stat_res.results_df.reset_index()
        logger.info(f"DESeq2 analysis completed. Found {len(results)} genes.")
        
        # Filter results based on significance
        significant_genes = results[
            (results['padj'] < par["padj_threshold"]) & 
            (results['log2FoldChange'].abs() > par["log2fc_threshold"])
        ]
        logger.info(f"Significant genes (padj < {par['padj_threshold']}, |log2FC| > {par['log2fc_threshold']}): {len(significant_genes)}")
        
        # Sort by log2FoldChange
        results = results.sort_values('log2FoldChange', ascending=False)
        
        # Add additional statistics
        results['abs_log2FoldChange'] = results['log2FoldChange'].abs()
        results['significant'] = (
            (results['padj'] < par["padj_threshold"]) & 
            (results['log2FoldChange'].abs() > par["log2fc_threshold"])
        )
    
    # Save results
    logger.info(f"Saving results to {par['output']}")
    results.to_csv(par["output"], index=False)
    
    # Log summary statistics
    upregulated = results[(results['log2FoldChange'] > 0) & results['significant']]
    downregulated = results[(results['log2FoldChange'] < 0) & results['significant']]
    
    logger.info(f"Summary:")
    logger.info(f"  Total genes analyzed: {len(results)}")
    logger.info(f"  Significant upregulated: {len(upregulated)}")
    logger.info(f"  Significant downregulated: {len(downregulated)}")

except Exception as e:
    logger.error(f"DESeq2 analysis failed: {str(e)}")
    logger.error(f"Analysis context: {contrast_column} contrast {contrast_values}, {len(mod.obs)} samples, {len(mod.var)} genes")
    logger.warning("Check your input data, design factors, and contrast specifications")
    
    raise e

logger.info("DESeq2 analysis completed successfully")