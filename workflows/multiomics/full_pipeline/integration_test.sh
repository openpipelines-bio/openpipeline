#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config

nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -entry test_wf2 \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config

nextflow \
  run . \
  -main-script workflows/multiomics/full_pipeline/main.nf \
  -entry test_wf3 \
  -resume \
  -profile docker,no_publish \
  -c workflows/utils/labels_ci.config

