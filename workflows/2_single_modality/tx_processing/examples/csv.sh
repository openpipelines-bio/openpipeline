#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

bin/nextflow \
  run . \
  -main-script workflows/2_single_modality/tx_processing/main.nf \
  --csv workflows/2_single_modality/tx_processing/examples/test_input.csv \
  --output output/ \
  --publishDir output_out/ \
  --logDir output_log/ \
  -resume \
  -c workflows/2_single_modality/tx_processing/nextflow.config