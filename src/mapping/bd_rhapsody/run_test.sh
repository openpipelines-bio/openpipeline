#!/bin/bash

echo ">> Checking whether requirement overrides work"
"./$meta_functionality_name" \
  --mode wta \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz"  \
  -r "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  -t "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf" \
  --putative_cell_call "mRNA" \
  --override_min_cores 1234 \
  --override_min_ram 5678 \
  --exact_cell_count 4900 \
  -o output2/ \
  --dryrun

if ! grep -q '"coresMin": 1234,' output2/pipeline.cwl; then
  echo Overriding minimum cores did not work
  exit 1
fi

if ! grep -q '"ramMin": 5678000,' output2/pipeline.cwl; then
  echo Overriding minimum ram did not work
  exit 1
fi


echo ">> Running $meta_functionality_name in WTA mode"
"$meta_executable" \
  --mode wta \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz"  \
  -r "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  -t "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf" \
  --putative_cell_call "mRNA" \
  --override_min_cores 1 \
  --override_min_ram 2 \
  --exact_cell_count 4900 \
  -o output/

echo ">> Checking whether output can be found"
[[ ! -f output/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1




echo ">> Running $meta_functionality_name in Targeted mode"
"$meta_executable" \
  --mode targeted \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-BCR_S1_L001_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-BCR_S1_L001_R2_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-mRNA_S5_L001_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-mRNA_S5_L001_R2_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-TCR_S3_L001_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_vdj/raw/RhapVDJDemo-TCR_S3_L001_R2_001_subset.fastq.gz" \
  --reference "$meta_resources_dir/bdrhap_vdj/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta" \
  --putative_cell_call "mRNA" \
  --override_min_cores 1 \
  --override_min_ram 2 \
  -o output_vdj/

echo ">> Checking whether output can be found"
[[ ! -f output_vdj/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1


echo ">>> Test finished successfully"
