#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

nextflow run . \
  -main-script workflows/multiomics/prot_multisample/main.nf \
  -profile docker,no_publish \
  -resume \
  -entry test_wf \
  -with-trace work/trace.txt

nextflow run . \
  -main-script workflows/multiomics/prot_multisample/main.nf \
  -profile docker,no_publish \
  -resume \
  -entry test2_wf \
  -with-trace work/trace.txt