#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# settings
ID=spaceranger_tiny_bcl
OUT="resources_test/$ID/"
DIR="$OUT"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# download bcl data
if [ ! -f "${OUT}/bcl/sample_sheet.csv" ]; then
  mkdir -p "$OUT/bcl"

  # download tar gz
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/supp/spatial-exp/spaceranger-tiny-bcl-1.0.0.tar.gz \
    --output "${OUT}/bcl/spaceranger-tiny-bcl-1.0.0.tar.gz"
  
  # untar
  tar -xf "${OUT}/bcl/spaceranger-tiny-bcl-1.0.0.tar.gz" \
    --strip-components=1 \
    -C "$OUT/bcl"

  # remove tar
  rm "${OUT}/bcl/spaceranger-tiny-bcl-1.0.0.tar.gz"

  # download sample sheet
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/supp/spatial-exp/spaceranger-tiny-bcl-simple-1.0.0.csv \
    --output "${OUT}/bcl/sample_sheet.csv"

fi

if [ ! -d "${OUT}/fastqs" ]; then
  mkdir -p "$OUT/fastqs"

  target/docker/demux/cellranger_mkfastq/cellranger_mkfastq \
    --input "${OUT}/bcl" \
    --sample_sheet "${OUT}/bcl/sample_sheet.csv" \
    --output "${OUT}/fastqs"
fi