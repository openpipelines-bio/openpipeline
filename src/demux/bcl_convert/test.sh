#!/bin/bash

# bcl-convert requires a v2 sample sheet
# bcl-convert is a bit more strict concerning filter files being present or not.
# We make a copy and make the necessary adaptations. 
# See workflows/resources_test_scripts/cellranger_tiny_bcl.sh for more information

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

echo ">>> Running executable"
$meta_executable \
  --input "$meta_resources_dir/bcl2" \
  --sample_sheet "$meta_resources_dir/bcl2/sample_sheet.csv" \
  --output "$TMPDIR/output1" \
  --test_mode true

echo ">>> Checking whether the output dir exists"
[[ ! -d "$TMPDIR/output1" ]] && echo "Output dir could not be found!" && exit 1

echo ">>> Checking whether output fastq files are created"
[[ ! -f  "$TMPDIR/output1/Undetermined_S0_R1_001.fastq.gz" ]] && echo "Output fastq files could not be found!" && exit 1
[[ ! -f "$TMPDIR/output1/s1_S1_R1_001.fastq.gz" ]] && echo "Output fastq files could not be found!" && exit 1

echo ">>> Test strict mode"
cp -r "$meta_resources_dir/bcl2" "$TMPDIR/strict_input"
rm "$TMPDIR/strict_input/Data/Intensities/L001/s_1_1101.clocs"

$meta_executable \
  --input "$TMPDIR/strict_input" \
  --sample_sheet "$meta_resources_dir/bcl2/sample_sheet.csv" \
  --output "$TMPDIR/output2" \
  --strict_mode false

echo ">>> Checking whether the output dir exists"
[[ ! -d "$TMPDIR/output2" ]] && echo "Output dir could not be found!" && exit 1

echo ">>> Checking whether output fastq files are created"
[[ ! -f  "$TMPDIR/output2/Undetermined_S0_R1_001.fastq.gz" ]] && echo "Output fastq files could not be found!" && exit 1
[[ ! -f "$TMPDIR/output2/s1_S1_R1_001.fastq.gz" ]] && echo "Output fastq files could not be found!" && exit 1

# print final message
echo ">>> Test finished successfully"

# do not remove this
# as otherwise your test might exit with a different exit code
exit 0
