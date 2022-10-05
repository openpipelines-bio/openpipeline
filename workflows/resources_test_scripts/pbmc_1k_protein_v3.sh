#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=pbmc_1k_protein_v3
OUT=resources_test/$ID/$ID
DIR=$(dirname "$OUT")

# ideally, this would be a versioned pipeline run
[ -d "$DIR" ] || mkdir -p "$DIR"

# dataset page:
# https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-gene-expression-and-cell-surface-protein-3-standard-3-0-0

# download metrics summary
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_metrics_summary.csv \
  --output "${OUT}_metrics_summary.csv"

# download counts h5 file
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  --output "${OUT}_filtered_feature_bc_matrix.h5"

# download counts matrix tar gz file
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.tar.gz \
  --output "${OUT}_filtered_feature_bc_matrix.tar.gz"

# extract matrix tar gz
mkdir -p "${OUT}_filtered_feature_bc_matrix"
tar -xvf "${OUT}_filtered_feature_bc_matrix.tar.gz" \
  -C "${OUT}_filtered_feature_bc_matrix" \
  --strip-components 1
rm "${OUT}_filtered_feature_bc_matrix.tar.gz"

# convert 10x h5 to h5mu
target/docker/convert/from_10xh5_to_h5mu/from_10xh5_to_h5mu \
  --input "${OUT}_filtered_feature_bc_matrix.h5" \
  --input_metrics_summary "${OUT}_metrics_summary.csv" \
  --output "${OUT}_filtered_feature_bc_matrix.h5mu" \
  --sample_id "$ID" \
  --id_to_obs_names true

# Convert h5mu to h5ad
python <<HEREDOC
import muon as mu
h5mu_data = mu.read_h5mu("${OUT}_filtered_feature_bc_matrix.h5mu")
h5mu_data.mod['rna'].write("${OUT}_filtered_feature_bc_matrix_rna.h5ad")
h5mu_data.mod['prot'].write("${OUT}_filtered_feature_bc_matrix_prot.h5ad")
HEREDOC

# run single sample
NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/multiomics/rna_singlesample/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_uss \
  --input "${OUT}_filtered_feature_bc_matrix.h5mu" \
  --output "`basename $OUT`_uss.h5mu" \
  --publishDir `dirname $OUT` \
  -resume

# run multisample
NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/multiomics/rna_multisample/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_ums \
  --input "${OUT}_uss.h5mu" \
  --output "`basename $OUT`_ums.h5mu" \
  --publishDir `dirname $OUT` \
  -resume

# run integration
NXF_VER=21.10.6 bin/nextflow \
  run . \
  -main-script workflows/multiomics/integration/main.nf \
  -profile docker \
  --id pbmc_1k_protein_v3_mms \
  --input "${OUT}_ums.h5mu" \
  --output "`basename $OUT`_mms.h5mu" \
  --publishDir `dirname $OUT` \
  --obs_covariates sample_id \
  -resume