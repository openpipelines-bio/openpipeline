#!/bin/bash

# collection of commands for testing workflows in repos

export NXF_VER=21.04.1

PAR_GENOME="http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA/GRCh38-PhiX-gencodev29/GRCh38-PhiX-gencodev29-20181205.tar.gz"
PAR_TRANSCRIPTOME="http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA/GRCh38-PhiX-gencodev29/gencodev29-20181205.gtf"
PAR_TSV="workflows/bd_rhapsody_wta/test/input.tsv"

nextflow drop https://github.com/openpipelines-bio/openpipeline.git

nextflow \
  run https://github.com/openpipelines-bio/openpipeline.git \
  -r 0.2.0 \
  -main-script workflows/1_ingestion/bd_rhapsody_wta/main.nf \
  -entry multi_wf \
  --tsv "$PAR_TSV" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  -resume
