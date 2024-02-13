#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

viash ns build -q rna_multisample

export NXF_VER=21.10.6

nextflow run . \
  -main-script src/workflows/rna/rna_multisample/test.nf \
  -profile docker,no_publish \
  -resume \
  -entry test_wf \
  -with-trace work/trace.txt \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config
