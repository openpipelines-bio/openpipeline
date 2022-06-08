#!/bin/bash

# collection of commands for testing workflows in repos

export NXF_VER=21.04.1
export VIASH_TEMP=$HOME/workspace/viash_temp/
export NXF_TEMP=$HOME/workspace/nxf_temp/

PAR_GENOME="resources_test/bd_rhapsody_wta_test/raw/GRCh38_primary_assembly_genome_chr1.tar.gz"
PAR_TRANSCRIPTOME="resources_test/bd_rhapsody_wta_test/raw/gencode_v40_annotation_chr1.gtf"
PAR_CSV='resources_test/bd_rhapsody_wta_test/raw/input.csv'

nextflow \
  run https://github.com/openpipelines-bio/openpipeline.git \
  -r 0.3.1 \
  -main-script workflows/1_ingestion/bd_rhapsody_wta/main.nf \
  --csv "$PAR_CSV" \
  --reference_genome "$PAR_GENOME" \
  --transcriptome_annotation "$PAR_TRANSCRIPTOME" \
  --output output/ \
  --publishDir output/ \
  -with-docker \
  -resume \
  -latest
