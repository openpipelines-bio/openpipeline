#!/bin/bash

# collection of commands for testing workflows in repos

export NXF_VER=21.04.1

PAR_GENOME="resources_test/bd_rhapsody_wta_test/raw/GRCh38-PhiX-gencodev29-20181205.tar.gz"
PAR_TRANSCRIPTOME="resources_test/bd_rhapsody_wta_test/raw/gencodev29-20181205.gtf"
PAR_CSV='resources_test/bd_rhapsody_wta_test/raw/input.csv'

nextflow drop https://github.com/openpipelines-bio/openpipeline.git

nextflow \
  run . \
  -main-script workflows/1_ingestion/bd_rhapsody_wta/main.nf \
  --csv "$PAR_CSV" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  -resume
