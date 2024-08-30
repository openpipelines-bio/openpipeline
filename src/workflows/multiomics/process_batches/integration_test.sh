#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

viash ns build -q process_batches

export NXF_VER=24.04.4

nextflow \
  run . \
  -main-script src/workflows/multiomics/process_batches/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  -resume


nextflow \
  run . \
  -main-script src/workflows/multiomics/process_batches/test.nf \
  -entry test_wf2 \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  -resume
