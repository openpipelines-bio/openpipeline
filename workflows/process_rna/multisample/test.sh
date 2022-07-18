#!/bin/bash

nextflow run . \
  -main-script workflows/process_rna/multisample/main.nf \
  -profile docker \
  --id foo \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_uss.h5mu \
  --publishDir output