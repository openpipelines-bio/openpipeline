#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

viash ns build -q multisample

export NXF_VER=23.04.2

nextflow \
  run . \
  -main-script src/workflows/multiomics/multisample/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -resume


nextflow \
  run . \
  -main-script src/workflows/multiomics/multisample/test.nf \
  -entry test_wf2 \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -resume
