#!/bin/bash

# collection of commands for testing workflows in repos

export NXF_VER=21.04.1

PAR_GENOME="resources_test/bd_rhapsody_wta_test/raw/GRCh38-PhiX-gencodev29-20181205.tar.gz"
PAR_TRANSCRIPTOME="resources_test/bd_rhapsody_wta_test/raw/gencodev29-20181205.gtf"
PAR_INPUTS='resources_test/bd_rhapsody_wta_test/raw/sample_R1_.fastq.gz;resources_test/bd_rhapsody_wta_test/raw/sample_R2_.fastq.gz'

nextflow drop https://github.com/openpipelines-bio/openpipeline.git

nextflow \
  run https://github.com/openpipelines-bio/openpipeline.git \
  -r main_build \
  -main-script workflows/1_ingestion/bd_rhapsody_wta/main.nf \
  --id "sample_test" \
  --input "$PAR_INPUTS" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  -resume \
  -latest
