#!/bin/bash

mkdir output

"./$meta_functionality_name" \
  -i "$meta_resources_dir/bd_rhapsody_wta_1_10_test/raw/sample_R1_.fastq.gz" \
  -i "$meta_resources_dir/bd_rhapsody_wta_1_10_test/raw/sample_R2_.fastq.gz"  \
  -r "$meta_resources_dir/bd_rhapsody_wta_test/raw/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  -t "$meta_resources_dir/bd_rhapsody_wta_test/raw/gencode_v40_annotation_chr1.gtf" \
  --subsample 0.2 \
  -o output/

[[ ! -f output/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1

echo ">>> Test finished successfully"
