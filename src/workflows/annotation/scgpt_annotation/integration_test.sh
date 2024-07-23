#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

viash ns build -q scgpt_annotation

nextflow run . \
  -main-script src/workflows/annotation/scgpt_annotation/test.nf \
  -profile docker,no_publish \
  -entry test_wf \
  -with-trace work/trace.txt \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config
