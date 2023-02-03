#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

nextflow \
  run . \
  -main-script workflows/ingestion/make_reference/main.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -resume \
  -with-trace work/trace.txt \
  -with-dag workflows/ingestion/make_reference/graph.dot
