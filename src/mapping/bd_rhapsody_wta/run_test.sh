#!/bin/bash

mkdir output

# temporarily decrease RAM and CPU requirements to be able to work on CI systems
sed -i 's#"ramMin": [^,]*,#"ramMin": 4000,#' "$meta_resources_dir/rhapsody_wta_1.9.1_nodocker.cwl"
sed -i 's#"coresMin": [^,]*,#"coresMin": 1,#' "$meta_resources_dir/rhapsody_wta_1.9.1_nodocker.cwl"

echo ">> Running $meta_functionality_name"
"./$meta_functionality_name" \
  -i "$meta_resources_dir/bd_rhapsody_wta_test/raw/sample_R1_.fastq.gz" \
  -i "$meta_resources_dir/bd_rhapsody_wta_test/raw/sample_R2_.fastq.gz"  \
  -r "$meta_resources_dir/bd_rhapsody_wta_test/raw/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  -t "$meta_resources_dir/bd_rhapsody_wta_test/raw/gencode_v40_annotation_chr1.gtf" \
  --subsample 0.2 \
  -o output/

echo ">> Checking whether output can be found"
[[ ! -f output/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1

echo ">>> Test finished successfully"
