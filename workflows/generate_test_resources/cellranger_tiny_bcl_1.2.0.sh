#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

OUT=resources_test/cellranger_tiny_bcl_1.2.0/
DIR="$OUT"
S3DIR=$(echo "$DIR" | sed 's#resources_test#s3://openpipelines-data#')

target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz \
  --output "${OUT}/cellranger-tiny-bcl-1.2.0.tar.gz"

target/docker/download/download_file/download_file \
  --input https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv \
  --output "${OUT}/cellranger-tiny-bcl-simple-1.2.0.csv"

# if subsampling is needed, uncomment the following:

# # extract tiny bcl tar
# tar -xf cellranger-tiny-bcl-1.2.0.tar.gz
# rm cellranger-tiny-bcl-1.2.0.tar.gz

# # keep only C1.1 to C3.1 basecalls
# mkdir cellranger-tiny-bcl-1.2.0/Data/Intensities/BaseCalls/L001new/
# cp -r cellranger-tiny-bcl-1.2.0/Data/Intensities/BaseCalls/L001/C[1-4].1 cellranger-tiny-bcl-1.2.0/Data/Intensities/BaseCalls/L001new/
# rm -r cellranger-tiny-bcl-1.2.0/Data/Intensities/BaseCalls/L001/
# mv cellranger-tiny-bcl-1.2.0/Data/Intensities/BaseCalls/L001new/ cellranger-tiny-bcl-1.2.0/Data/Intensities/BaseCalls/L001/

# # todo:
# # do I need to edit 'cellranger-tiny-bcl-1.2.0/Data/Intensities/config.xml' to reduce the number of cycles from 132 to 4?

# # create tar
# tar -czvf cellranger-tiny-bcl-1.2.0.tar.gz cellranger-tiny-bcl-1.2.0
# rm -r cellranger-tiny-bcl-1.2.0

aws s3 sync --profile di "$DIR" "$S3DIR"
