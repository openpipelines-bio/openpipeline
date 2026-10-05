#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

ID=geneformer
OUT=resources_test/$ID

# Geneformer-V1-10M is the smallest pretrained Geneformer model (~40 MB,
# hidden size 256, input size 2048). The tests use it to keep the download
# small; production runs use a V2 model (e.g. Geneformer-V2-316M).
# Pinned to the same commit as the geneformer package in
# src/base/geneformer_engine.yaml, so tokens and model match.
COMMIT=1f7fbae4e469a5f4f1af8c111a529cfe1b3829f5
MODEL_URL="https://huggingface.co/ctheodoris/Geneformer/resolve/$COMMIT/Geneformer-V1-10M"
MODEL_DIR="$OUT/Geneformer-V1-10M"

mkdir -p "$MODEL_DIR"
for file in config.json model.safetensors training_args.bin; do
  if [ ! -s "$MODEL_DIR/$file" ]; then
    echo "Downloading $file"
    wget -q "$MODEL_URL/$file" -O "$MODEL_DIR/$file"
  fi
done

# The workflow tests also use the shared pbmc_1k_protein_v3 dataset, which is
# created by resources_test_scripts/pbmc_1k_protein_v3.sh.

# aws s3 sync --profile di "$OUT" "s3://openpipelines-data/$ID" --delete --dryrun
