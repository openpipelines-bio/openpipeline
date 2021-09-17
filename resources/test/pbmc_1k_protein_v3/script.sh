#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUT=resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3

wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_metrics_summary.csv -O ${OUT}_metrics_summary.csv

target/docker/download/download_10x_dataset/download_10x_dataset \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5 \
  --output ${OUT}_raw_feature_bc_matrix.h5 \
  --min_library_size 1000 \
  --min_cells_per_gene 500 \
  --keep_feature_types "Antibody Capture"
  
target/docker/convert/convert_10x_to_h5ad/convert_10x_to_h5ad \
  --input ${OUT}_raw_feature_bc_matrix.h5 \
  --output ${OUT}_raw_feature_bc_matrix.h5ad \
  --compression gzip
