#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

bin/nextflow \
  run . \
  -main-script workflows/integration/mutliomics/main.nf \
  -entry test_wf \
  -resume \
  -profile docker \
  -c workflows/utils/labels_ci.config \
  --id "mouse;human" \
  --input "$REPO_ROOT/resources_test/concat/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset.h5mu;$REPO_ROOT/resources_test/concat/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu" \
  --publish_dir "foo/"
