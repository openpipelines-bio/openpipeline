#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.04.1

PAR_INPUT="resources_test/cellranger-tiny-bcl-1.2.0/cellranger-tiny-bcl-1.2.0.tar.gz"
PAR_SAMPLESHEET="resources_test/cellranger-tiny-bcl-1.2.0/cellranger-tiny-bcl-simple-1.2.0.csv"

nextflow run . \
  -main-script workflows/0_tmp/demultiplex-cellranger/demultiplex-cellranger.nf \
  --input "$PAR_INPUT" \
  --samplesheet "$PAR_SAMPLESHEET" \
  --publishDir output/ \
  -resume
