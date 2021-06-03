#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# collection of commands for testing workflows in repos

export NXF_VER=21.04.1

PAR_GENOME="http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA/GRCh38-PhiX-gencodev29/GRCh38-PhiX-gencodev29-20181205.tar.gz"
PAR_TRANSCRIPTOME="http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA/GRCh38-PhiX-gencodev29/gencodev29-20181205.gtf"
PAR_INPUTS="http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R1_.fastq.gz;http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R2_.fastq.gz"

nextflow drop https://github.com/openpipeline-bio/openpipeline.git

nextflow \
  run https://github.com/openpipeline-bio/openpipeline.git \
  -r main_build \
  -main-script workflows/bd_rhapsody_wta/main.nf \
  -entry single_wf \
  --id "sample_RSEC" \
  --input "$PAR_INPUTS" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  -resume
  
nextflow \
  run . \
  -main-script workflows/bd_rhapsody_wta/main.nf \
  -entry single_wf \
  --id "sample_RSEC" \
  --input "$PAR_INPUTS" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  -resume
