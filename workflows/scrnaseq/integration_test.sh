#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

nextflow run . \
  -main-script workflows/scrnaseq/main.nf \
  -profile docker,nopublish \
  -resume \
  -entry test_singlesample \
  -with-trace work/trace.txt

nextflow run . \
  -main-script workflows/scrnaseq/main.nf \
  -profile docker,nopublish \
  -resume \
  -entry test_multisample \
  -with-trace work/trace.txt