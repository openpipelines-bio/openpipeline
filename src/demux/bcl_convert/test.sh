#!/bin/bash

set -ex

# We are using the tiny bcl dataset provided by Illumina:
#   https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq
# Unfortunately,
#   1. the sample sheet delivered with it does not work with bcl-convert (v1 of the format)
#   2. 2 filter files are missing from the run directory that bcl-convert requires to run
#
# We worked around this by
#   1. Manually editing a sample sheet file suited for bcl-convert (format v2)
#   2. Adding a filter file
#
# The filter file is a binary file, we just created an empty file use that.
# bcl-convert might complain about it, but at least something is written out.
# An alternative is to use a filter file from a different project. This also generates
# a warning, but the fastq ouput files contain reads. The drawback is that those filter files
# are generally above 100MB in size.
#
# TODO: Check if a (binary) filter file can be generated that is small but works.

cp -r bcl local_copy
cp s_1_1101.filter local_copy/Data/Intensities/BaseCalls/L001/s_1_1101.filter

echo ">>> Running executable"
./bcl_convert \
  --input local_copy \
  --sample_sheet SampleSheet.csv \
  --output fastq \
  --reports reports \
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
