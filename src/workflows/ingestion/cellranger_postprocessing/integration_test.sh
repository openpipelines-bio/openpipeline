#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

nextflow \
  run . \
  -main-script src/workflows/ingestion/cellranger_postprocessing/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -with-trace work/trace.txt \
  -resume

nextflow \
  run . \
  -main-script src/workflows/ingestion/cellranger_postprocessing/test.nf \
  -entry test_wf2 \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -with-trace work/trace.txt \
  -resume
