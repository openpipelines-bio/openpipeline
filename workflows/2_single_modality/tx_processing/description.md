# Conversion

## Goal
Pipeline which takes h5ad files containing single cell tx-data, filters and performs a per sample processing.

![Per sample transcriptomics processing pipeline](../../../docs/figures/per-sample-tx.png)

_Overview of the per sample transcriptomics processing pipeline._

Input: 
- single sample h5ad object containing the raw counts
- standard filtering parameters
    - min_UMI_count
    - max_UMI_count
    - min_gene_count
    - max_gene_count
    - percent_mito
    - exclude_doublets
- optional filtering parameters

## Input

- ### ___input.h5ad*__

    - {raw}/X

- ### __leiden-resolutions__

List of the leiden resolutions to be performed on the dataset.


## Output

- ### {input_name}-converted.h5ad

    - X : lognorm reduced data
    - {raw}/X : raw expression counts

    - obsm/X_pca
    - obsm/X_umap
    
    - obs/leiden.res.{leiden-resolutions}
    The raw count data that was in the original files.
