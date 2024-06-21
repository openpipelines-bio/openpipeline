#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

viash ns build -q ingestion/cellranger_mapping --setup cb --platform nextflow

export NXF_VER=21.10.6

nextflow \
  run . \
  -main-script src/workflows/ingestion/cellranger_mapping/test.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  -with-trace work/trace.txt