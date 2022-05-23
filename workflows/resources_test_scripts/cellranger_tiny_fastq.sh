#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# settings
ID=cellranger_tiny_fastq
OUT="resources_test/$ID/"
DIR="$OUT"

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT


# download cellranger tar gz
cellranger_tar_gz="${OUT}/cellranger-6.1.2.tar.gz"
if [ ! -f "$cellranger_tar_gz" ]; then
  echo "Download Cell Ranger 6.1.2 manually first!"
  exit 1
fi

# untar fastqs
cellranger_tiny_fastq="${OUT}/cellranger_tiny_fastq"
if [ ! -f "${cellranger_tiny_fastq}/tinygex_S1_L001_R1_001.fastq.gz" ]; then
  mkdir -p "$cellranger_tiny_fastq"
  
  tar -xzf "$cellranger_tar_gz" \
    -C "$cellranger_tiny_fastq" \
    "cellranger-6.1.2/external/cellranger_tiny_fastq" \
    --strip-components=3
fi

# untar ref
cellranger_tiny_ref="${OUT}/cellranger_tiny_ref"
if [ ! -f "${cellranger_tiny_ref}/reference.json" ]; then
  mkdir -p "$cellranger_tiny_ref"
  
  tar -xzf "$cellranger_tar_gz" \
    -C "$cellranger_tiny_ref" \
    "cellranger-6.1.2/external/cellranger_tiny_ref" \
    --strip-components=3
fi

bam_dir="${OUT}/bam"
if [ ! -f "$bam_dir" ]; then
  mkdir -p "$bam_dir"

  target/docker/mapping/cellranger_count/cellranger_count \
    --input "$cellranger_tiny_fastq" \
    --reference "$cellranger_tiny_ref" \
    --output "$bam_dir"
fi

filtered_h5="${OUT}/filtered.h5"
if [ ! -f "$filtered_h5" ]; then
  target/docker/mapping/cellranger_count_split/cellranger_count_split \
    --input "$bam_dir" \
    --filtered_h5 "$filtered_h5"
fi

filtered_h5ad="${OUT}/filtered.h5ad"
if [ ! -f "$filtered_h5ad" ]; then
  target/docker/convert/from_10xh5_to_h5ad/from_10xh5_to_h5ad \
    --input "$filtered_h5" \
    --output "$filtered_h5ad"
fi

filtered_h5mu="${OUT}/filtered.h5mu"
if [ ! -f "$filtered_h5mu" ]; then
  target/docker/convert/from_10xh5_to_h5mu/from_10xh5_to_h5mu \
    --input "$filtered_h5" \
    --output "$filtered_h5mu"
fi
