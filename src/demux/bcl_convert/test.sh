#!/bin/bash



# bcl-convert requires a v2 sample sheet
# bcl-convert is a bit more strict concerning filter files being present or not.
# We make a copy and make the necessary adaptations. 
# See workflows/resources_test_scripts/cellranger_tiny_bcl.sh for more information

echo ">>> Running executable"
$meta_executable \
  --input "$meta_resources_dir/bcl2" \
  --sample_sheet "$meta_resources_dir/bcl2/sample_sheet.csv" \
  --output fastq \
  --test_mode true

echo ">>> Checking whether the output dir exists"
[[ ! -d fastq ]] && echo "Output dir could not be found!" && exit 1

echo ">>> Checking whether output fastq files are created"
[[ ! -f fastq/Undetermined_S0_R1_001.fastq.gz ]] && echo "Output fastq files could not be found!" && exit 1
[[ ! -f fastq/s1_S1_R1_001.fastq.gz ]] && echo "Output fastq files could not be found!" && exit 1

# print final message
echo ">>> Test finished successfully"

# do not remove this
# as otherwise your test might exit with a different exit code
exit 0
