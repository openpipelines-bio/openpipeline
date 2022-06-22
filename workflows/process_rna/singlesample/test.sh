#!/bin/bash

nextflow run . \
  -main-script workflows/2_unimodal_singlesample/rna/main.nf \
  -profile docker \
  -resume \
  --id foo \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu \
  --publishDir output