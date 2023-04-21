#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
OUT=resources_test/test_param_list.yaml

cat > $OUT << HERE
- id: "mouse"
  input: s3://openpipelines-data/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu
  publish_dir: "foo_remote/"
  obs_covariates: "sample_id"
  rna_min_counts: 2
  prot_min_counts: 3
- id: "human"
  input: "s3://openpipelines-data/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu
  publish_dir: "foo_remote/"
  obs_covariates: "sample_id"
  rna_min_counts: 2
  prot_min_counts: 3
HERE