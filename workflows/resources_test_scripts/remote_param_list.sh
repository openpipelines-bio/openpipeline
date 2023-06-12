#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
mkdir -p "resources_test/remote_param_list/"
OUT=resources_test/remote_param_list/test_param_list.yaml
OUT_CSV=resources_test/remote_param_list/test_param_list.csv
OUT_JSON=resources_test/remote_param_list/test_param_list.json

cat > $OUT << HERE
- id: "mouse"
  input: s3://openpipelines-data/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu
  publish_dir: "foo_remote/"
  rna_min_counts: 2
  prot_min_counts: 3
- id: "human"
  input: s3://openpipelines-data/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu
  publish_dir: "foo_remote/"
  rna_min_counts: 2
  prot_min_counts: 3
HERE

cat > $OUT_CSV << EOF
"id","input","publish_dir","rna_min_counts","prot_min_counts"
"mouse","s3://openpipelines-data/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu","foo_remote/","2","3"
"human","s3://openpipelines-data/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu","foo_remote/","2","3"
EOF

cat > $OUT_JSON << HERE
[
    {
        "id": "mouse",
        "input": "s3://openpipelines-data/concat_test_data/e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
        "publish_dir": "foo_remote/",
        "rna_min_counts": 2,
        "prot_min_counts": 3
    },
    {
        "id": "human",
        "input": "s3://openpipelines-data/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset_unique_obs.h5mu",
        "publish_dir": "foo_remote/",
        "rna_min_counts": 2,
        "prot_min_counts": 3
    }
]
HERE