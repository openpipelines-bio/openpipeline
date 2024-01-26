#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

viash ns build -q prot_multisample

export NXF_VER=21.10.6

nextflow run . \
  -main-script src/workflows/prot/prot_multisample/test.nf \
  -profile docker,no_publish \
  -entry test_wf \
  -with-trace work/trace.txt