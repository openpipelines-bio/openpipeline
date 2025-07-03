import scanpy as sc
import mudata as mu
import random
import numpy as np

## VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "mod": "rna",
    "min_cells": 100,
    "target_cluster_id": "donor_id",
    "grouping_column": "cell_type",  # Add this parameter
    "min_cells_per_sample": 50,     # Add this parameter
    "n_replicates": 1,              # Add this parameter
}

meta = {"resources_dir": "src/utils"}
### VIASH END

import sys
sys.path.append(meta["resources_dir"])

from setup_logger import setup_logger

logger = setup_logger()

mdata = mu.read_h5mu(par["input"])
mod = mdata.mod[par["mod"]]

sc.pp.filter_genes(
    mod,
    min_cells=par["min_cells"],
)

mod.obs_names_make_unique()
mod.var_names_make_unique()

mod.obs[par["target_cluster_id"]] = mod.obs[par["target_cluster_id"]].cat.remove_unused_categories()

mod_dge = mod.copy()

pbs = {}
for cluster in mod_dge.obs[par["target_cluster_id"]].cat.categories:
    logger.info(f"Processing cluster {cluster}")
    
    pbs[cluster] = []
    
    # Get cells from this cluster
    cluster_cells = mod_dge[mod_dge.obs[par["target_cluster_id"]] == cluster]
    
    grouping_values = cluster_cells.obs[par["grouping_column"]].unique()
    
    # Process each group within this cluster
    for group in grouping_values:
        group_cells = cluster_cells[cluster_cells.obs[par["grouping_column"]] == group]
        
        if len(group_cells) < par["min_cells_per_sample"]:
            logger.warning(f"Skipping {cluster}-{group}: only {len(group_cells)} cells (minimum: {par['min_cells_per_sample']})")
            continue
            
        # Create pseudoreplicates
        indices = list(group_cells.obs_names)
        random.shuffle(indices)
        indices_split = np.array_split(np.array(indices), par["n_replicates"])
        
        for i, pseudo_rep_indices in enumerate(indices_split):
            if len(pseudo_rep_indices) >= par["min_cells_per_sample"]:
                # Sum expression across cells to create pseudobulk sample
                pseudobulk_cells = group_cells[pseudo_rep_indices]
                
                rep_adata = sc.AnnData(
                    X=pseudobulk_cells.X.sum(axis=0),
                    var=pseudobulk_cells.var[[]]  # Empty var dataframe with same index
                )
                
                # Add metadata
                rep_adata.obs_names = [f'{group}_{cluster}_rep{i}']
                rep_adata.obs['replicate'] = i
                rep_adata.obs['group'] = str(group)
                rep_adata.obs['cluster'] = str(cluster)
                rep_adata.obs['n_cells'] = len(pseudo_rep_indices)
                
                # Add original metadata if available
                first_cell_metadata = pseudobulk_cells.obs.iloc[0]
                for col in pseudobulk_cells.obs.columns:
                    if col not in ['replicate', 'group', 'cluster', 'n_cells']:
                        rep_adata.obs[col] = str(first_cell_metadata[col])
                
                pbs[cluster].append(rep_adata)
                logger.info(f"Created pseudobulk sample: {group}_{cluster}_rep{i} with {len(pseudo_rep_indices)} cells")

# Concatenate all pseudobulk samples
logger.info("Concatenating pseudobulk samples...")
all_pseudobulk_samples = []
for cluster, samples in pbs.items():
    if samples:  # Only add if there are samples
        all_pseudobulk_samples.extend(samples)

if all_pseudobulk_samples:
    adata_pseudo = sc.concat(all_pseudobulk_samples)
    logger.info(f"Created {len(all_pseudobulk_samples)} pseudobulk samples with shape: {adata_pseudo.shape}")
    
    # Convert object columns to string for compatibility
    for col in adata_pseudo.obs.columns:
        if adata_pseudo.obs[col].dtype == 'object':
            adata_pseudo.obs[col] = adata_pseudo.obs[col].astype(str)
            
    logger.info("Pseudobulk generation completed successfully")
else:
    logger.error("No pseudobulk samples were generated!")
    
print("here")