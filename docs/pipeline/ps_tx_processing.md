# Per sample transcriptomics processing

The per sample pipeline takes raw count data from an h5ad object and filters based on the counts. 

## UMI, Gene, mito content and doublet filtering


![Per sample transcriptomics processing pipeline](../figures/per-sample-tx.png)

_Overview of the per sample transcriptomics processing pipeline._

Input: 
- (h5ad)/raw
- standard filtering parameters
    - min_UMI_count: int
    - max_UMI_count: int
    - min_gene_count: int
    - max_gene_count: int
    - percent_mito: float
    - exclude_doublets: boolean

Tools:
    - Scanpy 
    - Scrublet