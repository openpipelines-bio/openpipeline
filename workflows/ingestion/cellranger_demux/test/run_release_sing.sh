#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6
export NXF_SINGULARITY_CACHEDIR="$HOME/.cache/singularity"

bin/nextflow \
  run https://github.com/czbiohub/utilities \
  -r 1.0.0 \
  -main-script workflows/1_ingestion/cellranger_demux/main.nf \
  -resume \
  -latest \
  -with-singularity \
  --id tiny_bcl \
  --input resources_test/cellranger_tiny_bcl/bcl \
  --sample_sheet resources_test/cellranger_tiny_bcl/bcl/sample_sheet.csv \
  --publishDir temp