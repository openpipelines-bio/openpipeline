#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

viash ns build -q '^workflows'

nextflow \
  run . \
  -main-script src/workflows/multiomics/process_samples/test.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config

# Same as above but with remote yaml file.
nextflow \
  run . \
  -main-script src/workflows/multiomics/process_samples/test.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  --param_list s3://openpipelines-data/remote_param_list/test_param_list.yaml

# Same as above but with remote json file.
nextflow \
  run . \
  -main-script src/workflows/multiomics/process_samples/test.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  --param_list s3://openpipelines-data/remote_param_list/test_param_list.json

# Same as above but with remote csv file.
nextflow \
  run . \
  -main-script src/workflows/multiomics/process_samples/test.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  --param_list s3://openpipelines-data/remote_param_list/test_param_list.csv

nextflow \
  run . \
  -main-script src/workflows/multiomics/process_samples/test.nf \
  -entry test_wf2 \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config

nextflow \
  run . \
  -main-script src/workflows/multiomics/process_samples/test.nf \
  -entry test_wf3 \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config

nextflow \
  run . \
  -main-script src/workflows/multiomics/process_samples/test.nf \
  -entry test_wf4 \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config
