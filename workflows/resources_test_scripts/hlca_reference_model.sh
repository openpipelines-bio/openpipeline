#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=hlca_reference_model
OUT=resources_test/$ID/$ID
DIR=$(dirname "$OUT")

# ideally, this would be a versioned pipeline run
[ -d "$DIR" ] || mkdir -p "$DIR"

# download and unarchive pre-trained scANVI model
wget https://zenodo.org/record/6337966/files/HLCA_reference_model.zip \
  -O "${OUT}.zip"
unzip "${OUT}.zip"
rm "${OUT}.zip"

# Test query data
# Source publication: Delorey, Toni M., et al. “COVID-19 tissue atlases reveal SARS-CoV-2 pathology and cellular targets.” Nature 595.7865 (2021): 107-113.
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5230nnn/GSM5230027/suppl/GSM5230027_04-P103142-S149-R01_raw_feature_bc_matrix.h5.gz \
  -O "${OUT}_query_test.h5.gz"
gzip -d "${OUT}_query_test.h5.gz"

# convert 10x h5 to h5mu
target/docker/convert/from_h5ad_to_h5mu/from_h5ad_to_h5mu \
  --input "${OUT}_query_test.h5" \
  --output "${OUT}_query_test.h5mu"
