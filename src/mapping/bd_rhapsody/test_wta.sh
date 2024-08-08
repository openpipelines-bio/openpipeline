#!/bin/bash

set -e

echo ">> Running $meta_name in WTA mode"
"$meta_executable" \
  --mode wta \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz"  \
  -r "$meta_resources_dir/reference_gencodev41_chr1/reference_bd_rhapsody.tar.gz" \
  -t "$meta_resources_dir/reference_gencodev41_chr1/reference.gtf.gz" \
  --putative_cell_call "mRNA" \
  ---cpus 1 \
  ---memory 2gb \
  --exact_cell_count 4900 \
  -o output/

echo ">> Checking whether output can be found"
[[ ! -f output/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1

# todo: check whether tempdir is empty??


echo ">>> Test finished successfully"
