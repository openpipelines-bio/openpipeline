#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=22.10.3

nextflow \
  run . \
  -main-script src/workflows/ingestion/cellranger_multi/test.nf \
  -entry test_wf \
  -resume \
  -profile no_publish,docker \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  -with-trace work/trace.txt


nextflow \
  run . \
  -main-script src/workflows/ingestion/cellranger_multi/test.nf \
  -entry test_wf2 \
  -resume \
  -profile no_publish,docker \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  -with-trace work/trace.txt