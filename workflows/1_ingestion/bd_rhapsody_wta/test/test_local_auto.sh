#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"


export NXF_VER=21.04.1

PAR_GENOME="/scratch/gencodev38chr6.tar.gz"
PAR_TRANSCRIPTOME="/scratch/gencodev38chr6-filtered.gtf"
PAR_DIR="/app/project/CS000182"

nextflow \
  run . \
  -main-script workflows/1_ingestion/bd_rhapsody_wta/main.nf \
  -entry auto_wf \
  --dir "$PAR_DIR" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  -resume \
  -w /scratch/nxf_work_openpipeline
