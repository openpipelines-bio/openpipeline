#!/bin/bash

rsync -avz ndnd:/mnt/ibm_seq2/NextSeq-02/seqdata/211022_VH00444_61_AAAMYY2M5 ~/211022_VH00444_61_AAAMYY2M5

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# settings
ID=czb_test
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
    --input https://github.com/nf-core/test-datasets/raw/demultiplex/testdata/MiSeq/220422_M11111_0222_000000000-K9H97.tar.gz \
    --output "${OUT}/bcl/raw.tar.gz"
  
  # untar
  tar -xf "${OUT}/bcl/raw.tar.gz" \
    --strip-components=1 \
    -C "$OUT/bcl"

  # remove tar
  rm "${OUT}/bcl/raw.tar.gz"

  # download sample sheet
  target/docker/download/download_file/download_file \
    --input https://raw.githubusercontent.com/nf-core/test-datasets/demultiplex/testdata/MiSeq/SampleSheet.csv \
    --output "${OUT}/bcl/sample_sheet.csv"
fi

if [ ! -f "${OUT}/fastqs" ]; then
  mkdir -p "$OUT/fastqs"

  target/docker/demux/cellranger_mkfastq/cellranger_mkfastq \
    --input "${OUT}/bcl" \
    --sample_sheet "${OUT}/bcl/sample_sheet.csv" \
    --output "${OUT}/fastqs"
fi

# if [ ! -f "${OUT}/bam" ]; then
#   mkdir -p "$OUT/bam"

#   target/docker/mapping/cellranger_count/cellranger_count \
#     --input "${OUT}/fastqs" \
#     --reference "resources_test/reference/refdata-gex-GRCh38-2020-A" \
#     --chemistry "SC5P-PE" \
#     --output "${OUT}/bam"
# fi


