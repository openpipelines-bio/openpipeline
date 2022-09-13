#!/bin/bash

## VIASH START
meta_executable="viash run src/mapping/bd_rhapsody/config.vsh.yaml --"
meta_resources_dir="resources_test"
## VIASH END

echo ">> Checking whether requirement overrides work"
$meta_executable \
  --mode wta \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  -i "$meta_resources_dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz"  \
  -r "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/GRCh38_primary_assembly_genome_chr1.tar.gz" \
  -t "$meta_resources_dir/bdrhap_ref_gencodev40_chr1/gencode_v40_annotation_chr1.gtf" \
  --putative_cell_call "mRNA" \
  ---n_proc 11 \
  ---memory 56gb \
  --exact_cell_count 4900 \
  -o output2/ \
  --dryrun

if ! grep -q '"coresMin": 11,' output2/pipeline.cwl; then
  echo Overriding minimum cores did not work
  exit 1
fi

if ! grep -q '"ramMin": 57344,' output2/pipeline.cwl; then
  echo Overriding minimum ram did not work
  exit 1
fi


echo ">>> Test finished successfully"
