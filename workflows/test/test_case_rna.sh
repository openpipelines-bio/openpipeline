#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

# Standard 10x h5 conversion
nextflow \
  run . \
  -main-script workflows/1_ingestion/conversion/main.nf \
  --input_type 10x_h5 \
  --input ./resources_test/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  --output output/ \
  -resume \
  -c ./workflows/1_ingestion/conversion/nextflow.config


# Run per sample tx_processing pipeline
nextflow \
  run . \
  -main-script workflows/2_single_modality/tx_processing/main.nf \
  --input ./output/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5.h5ad \
  --output output/ \
  -resume \
  -c ./workflows/2_single_modality/tx_processing/nextflow.config
