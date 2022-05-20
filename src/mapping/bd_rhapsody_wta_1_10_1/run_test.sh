#!/bin/bash

mkdir output

"./$meta_functionality_name" \
  -i "$meta_resources_dir/bd_rhapsody_wta_test/raw/sample_R1_.fastq.gz" \
  -i "$meta_resources_dir/bd_rhapsody_wta_test/raw/sample_R2_.fastq.gz"  \
  -r "$meta_resources_dir/bd_rhapsody_wta_test/raw/GRCh38-PhiX-gencodev29-20181205.tar.gz" \
  -t "$meta_resources_dir/bd_rhapsody_wta_test/raw/gencodev29-20181205.gtf" \
  --subsample 0.2 \
  -o output/

[[ ! -f output/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1

echo ">>> Test finished successfully"
