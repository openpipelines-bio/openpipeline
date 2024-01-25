#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=23.04.2

nextflow run . \
  -main-script src/workflows/multiomics/integration/totalvi_leiden/test.nf \
  -profile docker,no_publish \
  -entry test_wf \
  -c src/workflows/utils/labels_ci.config
