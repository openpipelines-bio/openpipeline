#!/bin/bash



# We are using the tiny bcl dataset provided by Illumina:
#   https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq
# Unfortunately,
#   1. the sample sheet delivered with it does not work with bcl2fastq (v1 of the format)
#   2. 2 filter files are missing from the run directory that bcl2fastq requires to run
#
# We worked around this by ignoring all missing entries

echo ">>> Running executable"
./bcl2fastq \
  --input bcl \
  --sample_sheet SampleSheet.csv \
  --output fastq \
  --separate_reports true \
  --reports reports \
  --ignore_missing

echo ">>> Checking whether the output dir exists"
[[ ! -d fastq ]] && echo "Output dir could not be found!" && exit 1

echo ">>> Checking whether output fastq files are created"
[[ ! -f fastq/Undetermined_S0_L001_R1_001.fastq.gz ]] && echo "Output fastq files could not be found!" && exit 1

echo ">>> Checking whether reports files are created"
[[ ! -d reports/html ]] && echo "Output reports files could not be found!" && exit 1

echo ">>> Running executable without separate reports"
./bcl2fastq \
  --input bcl \
  --sample_sheet SampleSheet.csv \
  --output fastq1 \
  --ignore_missing

echo ">>> Checking whether the output dir exists"
[[ ! -d fastq1 ]] && echo "Output dir could not be found!" && exit 1

echo ">>> Checking whether output fastq files are created"
[[ ! -f fastq1/Undetermined_S0_L001_R1_001.fastq.gz ]] && echo "Output fastq files could not be found!" && exit 1

echo ">>> Checking whether reports files are created"
[[ ! -d fastq1/Reports/html ]] && echo "Output reports files could not be found!" && exit 1

# print final message
echo ">>> Test finished successfully"

# do not remove this
# as otherwise your test might exit with a different exit code
exit 0
