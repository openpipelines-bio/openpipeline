#!/bin/bash

set -e

echo ">> Running $meta_name in Targeted mode"
"$meta_executable" \
  --mode targeted \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-BCR_S1_L001_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-BCR_S1_L001_R2_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-mRNA_S5_L001_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-mRNA_S5_L001_R2_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-TCR_S3_L001_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-TCR_S3_L001_R2_001_subset.fastq.gz" \
  --reference "$meta_resources_dir/bdrhap_vdj/raw/BD_Rhapsody_Immune_Response_Panel_Hs.fasta" \
  --putative_cell_call "mRNA" \
  ---cpus 1 \
  ---memory 2gb \
  -o output_vdj/

echo ">> Checking whether output can be found"
[[ ! -f output_vdj/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1

# todo: check whether tempdir is empty??

echo ">>> Test finished successfully"
