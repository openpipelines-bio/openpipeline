#!/bin/bash

set -eo pipefail

# Get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# Ensure the script runs from the root of the repository
cd "$REPO_ROOT"

# Settings
ID=spaceranger_tiny_fastq
OUT="resources_test/$ID/"

# Create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")

function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# Download fastq data
if [[ ! -f "${OUT}/fastq/sample_sheet.csv" ]]; then
  mkdir -p "$OUT/fastq"

  # Download fastqs
  target/docker/download/download_file/download_file \
    --input https://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_fastqs.tar \
    --output "${OUT}/fastq/V1_Adult_Mouse_Brain_fastqs.tar"

  # Untar
  tar -xf "${OUT}/fastq/V1_Adult_Mouse_Brain_fastqs.tar" \
    --strip-components=1 \
    -C "$OUT/fastq"

  # Remove tar
  rm "${OUT}/fastq/V1_Adult_Mouse_Brain_fastqs.tar"
fi

# Download image data
if [[ ! -f "${OUT}/image/V1_Adult_Mouse_Brain_image.tif" ]]; then
  mkdir -p "$OUT/image"

  # Download image
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Adult_Mouse_Brain/V1_Adult_Mouse_Brain_image.tif \
    --output "${OUT}/image/V1_Adult_Mouse_Brain_image.tif"
fi

# Download reference data
if [[ ! -f "${OUT}/reference/refdata-gex-mm10-2020-A" ]]; then
  mkdir -p "$OUT/reference"

  # Download reference
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/supp/spatial-exp/refdata-gex-mm10-2020-A.tar.gz \
    --output "${OUT}/reference/refdata-gex-mm10-2020-A.tar.gz"

  # Untar
  tar -xf "${OUT}/reference/refdata-gex-mm10-2020-A.tar.gz" \
    --strip-components=1 \
    -C "$OUT/reference"

  # Remove tar
  rm "${OUT}/reference/refdata-gex-mm10-2020-A.tar.gz"
fi

# Run spaceranger_count
if [[ ! -d "${OUT}/fastqs" ]]; then
  mkdir -p "$OUT/fastqs"

  target/docker/mapping/spaceranger_count/spaceranger_count \
    --input="${OUT}/fastqs" \
    --reference="${OUT}/reference/refdata-gex-mm10-2020-A" \
    --image="${OUT}/image/V1_Adult_Mouse_Brain_image.tif" \
    --slide="V19L01-041" \
    --area="C1"
fi
