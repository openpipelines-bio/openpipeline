#!/bin/bash

set -eo pipefail


# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=demuxafy_test_data
OUT=resources_test/$ID

# https://www.dropbox.com/s/m8u61jn4i1mcktp/TestData4PipelineSmall.tar.gz


wget "https://www.dropbox.com/s/m8u61jn4i1mcktp/TestData4PipelineSmall.tar.gz" -O ".cache/TestData4PipelineSmall.tar.gz"
tar -xvzf .cache/TestData4PipelineSmall.tar.gz --directory "./.cache/"
mkdir -p ./resources_test/demuxafy_test_data/

mv -vf .cache/TestData4PipelineSmall/test_dataset.vcf ./resources_test/demuxafy_test_data/
mv -vf .cache/TestData4PipelineSmall/test_dataset/outs/pooled.sorted.bam ./resources_test/demuxafy_test_data/