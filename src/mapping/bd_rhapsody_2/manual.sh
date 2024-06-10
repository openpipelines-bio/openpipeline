#!/bin/bash

dir=resources_test

viash run src/mapping/bd_rhapsody_2/config.vsh.yaml -- \
  --reads "$dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R1_001_subset.fastq.gz" \
  --reads "$dir/bdrhap_5kjrt/raw/12WTA_S1_L432_R2_001_subset.fastq.gz" \
  --reference_archive "/home/rcannood/Downloads/bdrhap/RhapRef_Human_WTA_2023-02.tar.gz" \
  --output output/bdrhap_5kjrt \
  ---cpus 30 \
  ---memory 60gb


  # --dryrun \