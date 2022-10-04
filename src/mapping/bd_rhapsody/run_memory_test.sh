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
  -r "$meta_resources_dir/reference_gencodev41_chr1/reference_bd_rhapsody.tar.gz" \
  -t "$meta_resources_dir/reference_gencodev41_chr1/reference.gtf.gz" \
  --putative_cell_call "mRNA" \
  ---cpus 11 \
  ---memory 56gb \
  --exact_cell_count 4900 \
  -o output2/ \
  --dryrun

if ! grep -q '"coresMin": 11,' output2/pipeline.cwl; then
  echo Overriding minimum cores did not work
  exit 1
fi
# (56 - 2) * 1024 = 55344
if ! grep -q '"ramMin": 55344,' output2/pipeline.cwl; then
  echo Overriding minimum ram did not work
  exit 1
fi


echo ">>> Test finished successfully"
