#!/bin/bash

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
# We make a copy and make the necessary adaptations. Please refer to the
# test included in bcl_convert itself for more information.
if [ ! -f "${OUT}/bcl2/RunParameters.xml" ]; then
  cp ${OUT}/bcl2/sample_sheet.csv ${OUT}/bcl2/sample_sheet.csv.backup
  cp -r ${OUT}/bcl/* "${OUT}/bcl2/"
  mv ${OUT}/bcl2/sample_sheet.csv.backup ${OUT}/bcl2/sample_sheet.csv
fi
if [ ! -f "${OUT}/bcl2/Data/Intensities/BaseCalls/L001/s_1_1101.filter" ]; then
  touch "${OUT}/bcl2/Data/Intensities/BaseCalls/L001/s_1_1101.filter"
fi

# if [ ! -f "${OUT}/bam" ]; then
#   mkdir -p "$OUT/bam"

#   target/docker/mapping/cellranger_count/cellranger_count \
#     --input "${OUT}/fastqs" \
#     --reference "resources_test/reference/refdata-gex-GRCh38-2020-A" \
#     --chemistry "SC5P-PE" \
#     --output "${OUT}/bam"
# fi

