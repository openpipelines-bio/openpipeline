#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# settings
ID=cellranger_tiny_bcl
OUT="resources_test/$ID/"
DIR="$OUT"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# download bcl data
if [ ! -f "${OUT}/bcl/sample_sheet.csv" ]; then
  mkdir -p "$OUT/bcl"

  # download tar gz
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz \
    --output "${OUT}/bcl/cellranger-tiny-bcl-1.2.0.tar.gz"
  
  # untar
  tar -xf "${OUT}/bcl/cellranger-tiny-bcl-1.2.0.tar.gz" \
    --strip-components=1 \
    -C "$OUT/bcl"

  # remove tar
  rm "${OUT}/bcl/cellranger-tiny-bcl-1.2.0.tar.gz"

  # download sample sheet
  target/docker/download/download_file/download_file \
    --input https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv \
    --output "${OUT}/bcl/sample_sheet.csv"
fi

if [ ! -d "${OUT}/fastqs" ]; then
  mkdir -p "$OUT/fastqs"

  target/docker/demux/cellranger_mkfastq/cellranger_mkfastq \
    --input "${OUT}/bcl" \
    --sample_sheet "${OUT}/bcl/sample_sheet.csv" \
    --output "${OUT}/fastqs"
fi

# bcl-convert requires a v2 sample sheet
# bcl-convert is a bit more strict concerning filter files being present or not.
# We make a copy and make the necessary adaptations. 

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

if [ ! -f "${OUT}/bcl2/sample_sheet.csv" ]; then
  mkdir "${OUT}/bcl2/"
  cp -r ${OUT}/bcl/* "${OUT}/bcl2/"
  cat > "${OUT}/bcl2/sample_sheet.csv" << HERE
[Header],,,,,,,,,
FileFormatVersion,2,,,,,,
RunName,hiseq_test,,,,,,
InstrumentPlatform,NextSeq,,,,,,
IndexOrientation,Forward,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
Read1Cycles,26,,,,,,,,,
Read2Cycles,98,,,,,,
,,,,,,,,,
[Sequencing_Settings],,,,,,,
,,,,,,,
[BCLConvert_Settings],,,,,,,
SoftwareVersion,3.8.4,,,,,,
NoLaneSplitting,true,,,,,,
FastqCompressionFormat,gzip,,,,,,
,,,,,,,,,
[BCLConvert_Data],,,,,,,
Sample_ID,index,,,,,,
s1,GGTTTACT,,,,,,
,,,,,,,
[Cloud_Settings],,,,,,,
GeneratedVersion,1.3.0.202111171923,,,,,,
,,,,,,,
[Cloud_Data],,,,,,,
Sample_ID,ProjectName,LibraryName,LibraryPrepKitName,IndexAdapterKitName,I7_Index_ID,Sample_Name,Description,Instrument,Type
s1,p1,s1_SI-P03-C9,,,IDT01,SI-P03-C9,s1,NextSeq,HighOutput_75cycles
HERE
  
  touch "${OUT}/bcl2/Data/Intensities/BaseCalls/L001/s_1_1101.filter"
fi
