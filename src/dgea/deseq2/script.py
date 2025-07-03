import scanpy as sc
import mudata as mu
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

## VIASH START
par = {
    "input": "pseudobulk_samples.h5mu",
    "output": "deseq2_results.csv",
    "design_factors": ["disease", "treatment"],
    "contrast_column": "treatment",
    "contrast_baseline": "ctrl",
    "contrast_treatment": "stim",
    "filter_genes_min_samples": None,
    "padj_threshold": 0.05,
    "log2fc_threshold": 0.0,
    "filter_gene_patterns": ["MIR\\d+", "AL\\d+", "LINC\\d+", "AC\\d+", "AP\\d+"],
    "var_gene_name" : "feature_name",
}

meta = {"resources_dir": "src/utils"}
### VIASH END

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

# Validate required columns exist
required_columns = par["design_factors"] + [par["contrast_column"]]
missing_columns = [col for col in required_columns if col not in mod.obs.columns]
if missing_columns:
    raise ValueError(f"Missing required columns in metadata: {missing_columns}")

# Prepare counts matrix
logger.info("Preparing counts matrix for DESeq2")
counts = pd.DataFrame(
    mod.X, 
    columns=mod.var_names,
    index=mod.obs_names
)

# Ensure counts are integers (required for DESeq2)
counts = counts.astype(int)

metadata = mod.obs.copy()

logger.info(f"Counts matrix shape: {counts.shape}")
logger.info(f"Design factors: {par['design_factors']}")

# Check contrast values exist
contrast_values = metadata[par["contrast_column"]].unique()
logger.info(f"Available values in {par['contrast_column']}: {contrast_values}")

if par["contrast_baseline"] not in contrast_values:
    raise ValueError(f"Baseline '{par['contrast_baseline']}' not found in {par['contrast_column']}")
if par["contrast_treatment"] not in contrast_values:
    raise ValueError(f"Treatment '{par['contrast_treatment']}' not found in {par['contrast_column']}")

# Create DESeq2 dataset
logger.info("Creating DESeq2 dataset")
try:
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors=par["design_factors"],
    )
    
    # Optional gene filtering
    if par["filter_genes_min_samples"] is not None:
        logger.info(f"Filtering genes expressed in at least {par['filter_genes_min_samples']} samples")
        sc.pp.filter_genes(dds, min_cells=par["filter_genes_min_samples"])
    else:
        # Default: filter genes expressed in all samples
        sc.pp.filter_genes(dds, min_cells=len(mod.obs))
    
    logger.info(f"After filtering: {dds.shape}")
    
    # Run DESeq2 analysis
    logger.info("Running DESeq2 analysis...")
    dds.deseq2()
    
    # Perform statistical test
    logger.info(f"Performing contrast: {par['contrast_column']} {par['contrast_treatment']} vs {par['contrast_baseline']}")
    stat_res = DeseqStats(
        dds, 
        contrast=(par["contrast_column"], par["contrast_treatment"], par["contrast_baseline"])
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
    
    # Filter out unwanted gene patterns if specified
    if par["filter_gene_patterns"]:
        pattern_string = "|".join(par["filter_gene_patterns"])
        before_filter = len(results)
        results_filtered = results[~results[par["var_gene_name"]].str.contains(pattern_string, na=False)]
        logger.info(f"Filtered out genes matching patterns {par['filter_gene_patterns']}: {before_filter} -> {len(results_filtered)}")
        results = results_filtered
    
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
    logger.info(f"  Top upregulated gene: {upregulated.iloc[0][par['var_gene_name']] if len(upregulated) > 0 else 'None'}")
    logger.info(f"  Top downregulated gene: {downregulated.iloc[-1][par['var_gene_name']] if len(downregulated) > 0 else 'None'}")

except Exception as e:
    logger.error(f"DESeq2 analysis failed: {str(e)}")
    logger.error(f"Analysis context: {par['contrast_column']} contrast '{par['contrast_treatment']}' vs '{par['contrast_baseline']}', {len(mod.obs)} samples, {len(mod.var)} genes")
    logger.warning("Check your input data, design factors, and contrast specifications")
    
    raise e

logger.info("DESeq2 analysis completed successfully")