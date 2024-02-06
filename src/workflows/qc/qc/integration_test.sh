#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

viash ns build -q '^workflows/qc/qc'

nextflow \
  run . \
  -main-script src/workflows/qc/qc/test.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config
