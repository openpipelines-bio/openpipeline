#!/bin/bash

nextflow run . \
  -main-script workflows/integration/multimodal_integration/main.nf \
  -profile docker \
  --id foo \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_ums.h5mu \
  --publishDir output