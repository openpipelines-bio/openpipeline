#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUT=resources/test/pbmc_1k_protein_v3/pbmc_1k_protein_v3
DIR=`dirname "$OUT"`
S3DIR=`echo "$DIR" | sed 's#resources/test#s3://openpipelines-data#'`

# ideally, this would be a versioned pipeline run

target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_metrics_summary.csv \
  --output "${OUT}_metrics_summary.csv"

target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  --output "${OUT}_filtered_feature_bc_matrix.h5"

target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tar.gz \
  --output "${OUT}_filtered_feature_bc_matrix.tar.gz"

mkdir -p "${OUT}_filtered_feature_bc_matrix"
tar -xvf "${OUT}_filtered_feature_bc_matrix.tar.gz" \
  -C "${OUT}_filtered_feature_bc_matrix" \
  --strip-components 1
rm "${OUT}_filtered_feature_bc_matrix.tar.gz"
  
target/docker/convert/convert_10x_h5_to_h5mu/convert_10x_h5_to_h5mu \
  --input "${OUT}_filtered_feature_bc_matrix.h5" \
  --output "${OUT}_filtered_feature_bc_matrix.h5mu"

aws s3 sync --profile xxx "$DIR" "$S3DIR"