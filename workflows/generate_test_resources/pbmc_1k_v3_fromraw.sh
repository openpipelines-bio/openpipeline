#!/bin/bash

# settings
ID=pbmc_1k_v3_fromraw
OUT=resources/test/$ID/$ID
DIR=$(dirname "$OUT")
S3DIR=$(echo "$DIR" | sed 's#resources/test#s3://openpipelines-data#')

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# create tempdir
TMPDIR=$(mktemp -d "$VIASH_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# download fastqs
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar \
  --output "$TMPDIR/pbmc_1k_v3_fastqs.tar"
tar -xf "$TMPDIR/pbmc_1k_v3_fastqs.tar" -C "$TMPDIR"

# download reference
target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz \
  --output "$TMPDIR/refdata-cellranger-GRCh38-3.0.0.tar.gz"
tar -xf "$TMPDIR/refdata-cellranger-GRCh38-3.0.0.tar.gz" -C "$TMPDIR"

# map
target/docker/mapping/cellranger_count/cellranger_count \
  --input "$TMPDIR/pbmc_1k_v3_fastqs/" \
  --transcriptome "$TMPDIR/refdata-cellranger-GRCh38-3.0.0/" \
  --output "$TMPDIR/output"

# convert
target/docker/convert/convert_10x_h5_to_h5mu/convert_10x_h5_to_h5mu \
  --input "$TMPDIR/output/filtered_feature_bc_matrix.h5" \
  --output "$TMPDIR/pbmc_1k_v3.h5mu"

# aws s3 sync --profile xxx "$DIR" "$S3DIR"
