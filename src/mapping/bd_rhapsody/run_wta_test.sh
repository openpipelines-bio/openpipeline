#!/bin/bash



echo ">> Running $meta_functionality_name in WTA mode"
"$meta_executable" \
  --mode wta \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz"  \
  -r "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  -t "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf" \
  --putative_cell_call "mRNA" \
  ---n_proc 1 \
  ---memory 2gb \
  --exact_cell_count 4900 \
  -o output/

echo ">> Checking whether output can be found"
[[ ! -f output/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1

# todo: check whether tempdir is empty??


echo ">>> Test finished successfully"
