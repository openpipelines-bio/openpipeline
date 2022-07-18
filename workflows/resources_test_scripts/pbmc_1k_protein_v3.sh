#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUT=resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3
DIR=$(dirname "$OUT")
S3DIR=$(echo "$DIR" | sed 's#resources_test#s3://openpipelines-data#')

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

target/docker/convert/from_10xh5_to_h5mu/from_10xh5_to_h5mu \
  --input "${OUT}_filtered_feature_bc_matrix.h5" \
  --output "${OUT}_filtered_feature_bc_matrix.h5mu"

target/docker/convert/from_10xh5_to_h5ad/from_10xh5_to_h5ad \
  --input "${OUT}_filtered_feature_bc_matrix.h5" \
  --output "${OUT}_filtered_feature_bc_matrix.h5ad"

NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/process_rna/singlesample/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_uss \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu \
  --publishDir resources_test/pbmc_1k_protein_v3/ \
  -resume

NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/process_rna/multisample/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_ums \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_uss.h5mu \
  --publishDir resources_test/pbmc_1k_protein_v3/ \
  -resume

NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/integration/multimodal_integration/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_mms \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_ums.h5mu \
  --publishDir resources_test/pbmc_1k_protein_v3/ \
  -resume