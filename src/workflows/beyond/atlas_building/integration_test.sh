#!/bin/bash
set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# Generate test data if needed
if [[ ! -f resources_test/beyond_test_data/donor_01.h5mu ]]; then
  bash resources_test_scripts/beyond_atlas_building_test_data.sh
fi

nextflow run . \
  -main-script src/workflows/beyond/atlas_building/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config
