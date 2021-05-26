#!/bin/bash

export NXF_VER=20.12.1-edge

PAR_GENOME="http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA/GRCh38-PhiX-gencodev29/GRCh38-PhiX-gencodev29-20181205.tar.gz"
PAR_TRANSCRIPTOME="http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA/GRCh38-PhiX-gencodev29/gencodev29-20181205.gtf"
PAR_INPUTS="http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R1_.fastq.gz;http://bd-rhapsody-public.s3.amazonaws.com/Rhapsody-WTA-test-data/sample_R2_.fastq.gz"

nextflow drop https://github.com/openpipeline-bio/openpipeline.git

nextflow \
  run https://github.com/openpipeline-bio/openpipeline.git \
  -r target_main \
  -latest \
  -main-script workflows/bd_rhapsody_wta/main.nf \
  -entry bd_rhapsody_wta_wf \
  --id "sample_RSEC" \
  --input "$PAR_INPUTS" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  -resume
  
  
  
  
nextflow \
  run https://github.com/openpipeline-bio/openpipeline.git \
  -r main \
  -latest \
  -main-script workflows/debug.nf
  
