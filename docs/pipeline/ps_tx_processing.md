# Per sample transcriptomics processing

The per sample pipeline takes raw count data from an h5ad object and filters based on the counts. 


![Per sample transcriptomics processing pipeline](../figures/per-sample-tx.png)
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

## UMI, Gene and mito content filtering

Tooling: Scanpy

## Doublet filtering

Tooling: Scrublet