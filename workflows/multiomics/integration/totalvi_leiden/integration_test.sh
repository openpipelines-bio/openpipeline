#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=23.04.2

nextflow run . \
  -main-script workflows/multiomics/integration/totalvi_leiden/main.nf \
  -profile docker,no_publish \
  -entry test_wf \
  -with-trace work/trace.txt \
  -with-dag workflows/multiomics/integration/totalvi_leiden/graph.dot \
  -c workflows/utils/labels_ci.config
