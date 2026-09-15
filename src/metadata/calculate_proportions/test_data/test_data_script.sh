#!/usr/bin/env bash
# Test data for metadata/calculate_proportions.
# Delegates to resources_test_scripts/beyond_trajectory_test_data.sh.
# Input: resources_test/beyond_test_data/atlas.h5mu
#   (has obs["participant_id"] and obs["subpopulation"])
#
# Run from the repo root:
#   bash src/metadata/calculate_proportions/test_data/test_data_script.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
BEYOND_DIR="${REPO_ROOT}/resources_test/beyond_test_data"

if [[ ! -f "${BEYOND_DIR}/atlas.h5mu" ]]; then
  echo "Generating BEYOND test data..."
  bash "${REPO_ROOT}/resources_test_scripts/beyond_trajectory_test_data.sh"
fi

echo "Test data for metadata/calculate_proportions:"
echo "  input:  ${BEYOND_DIR}/atlas.h5mu"
echo ""
echo "Run the component manually:"
echo ""
echo "  viash run src/metadata/calculate_proportions/config.vsh.yaml -- \\"
echo "    --input  ${BEYOND_DIR}/atlas.h5mu \\"
echo "    --output ${BEYOND_DIR}/proportions_output.h5mu"
