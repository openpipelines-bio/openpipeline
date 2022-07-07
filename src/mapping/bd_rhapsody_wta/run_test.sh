#!/bin/bash

echo ">> Checking whether requirement overrides work"
"./$meta_functionality_name" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz"  \
  -r "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  -t "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf" \
  --abseq_reference "$meta_resources_dir/bdrhap_5kjrt/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta" \
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


echo ">> Running $meta_functionality_name"
"$meta_executable" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz"  \
  -r "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  -t "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf" \
  --abseq_reference "$meta_resources_dir/bdrhap_5kjrt/raw/BDAbSeq_ImmuneDiscoveryPanel.fasta" \
  --putative_cell_call "mRNA" \
  --override_min_cores 1 \
  --override_min_ram 2 \
  --exact_cell_count 4900 \
  -o output/

echo ">> Checking whether output can be found"
[[ ! -f output/sample_RSEC_ReadsPerCell_Unfiltered.csv.gz ]] && echo "Output file could not be found!" && exit 1


echo ">>> Test finished successfully"
