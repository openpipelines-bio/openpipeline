#!/bin/bash

# collection of commands for testing workflows in repos

export NXF_VER=21.04.1

PAR_GENOME="resources_test/bd_rhapsody_wta_test/raw/GRCh38_primary_assembly_genome_chr1.tar.gz"
PAR_TRANSCRIPTOME="resources_test/bd_rhapsody_wta_test/raw/gencode_v40_annotation_chr1.gtf"
PAR_INPUTS='resources_test/bd_rhapsody_wta_test/raw/sample_R1_.fastq.gz;resources_test/bd_rhapsody_wta_test/raw/sample_R2_.fastq.gz'

nextflow \
  run . \
  -main-script workflows/ingestion/bd_rhapsody_wta/main.nf \
  --id "sample_test" \
  --input "$PAR_INPUTS" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  --publishDir output/ \
  -with-docker \
  -resume \
  -latest
