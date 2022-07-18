#!/bin/bash

nextflow run . \
  -main-script workflows/process_rna/singlesample/main.nf \
  -profile docker \
  -resume \
  --id foo \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu \
  --publishDir output