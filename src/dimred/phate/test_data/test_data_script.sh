#!/usr/bin/env bash
# Test data for dimred/phate.
# Delegates to resources_test_scripts/beyond_trajectory_test_data.sh.
# Input: resources_test/beyond_test_data/proportions_output.h5mu
#   (has obsm["proportions"] - participant-proportions broadcast per cell,
#    as produced by metadata/calculate_proportions)
#
# Run from the repo root:
#   bash src/dimred/phate/test_data/test_data_script.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
BEYOND_DIR="${REPO_ROOT}/resources_test/beyond_test_data"

if [[ ! -f "${BEYOND_DIR}/proportions_output.h5mu" ]]; then
  echo "Generating BEYOND test data..."
  bash "${REPO_ROOT}/resources_test_scripts/beyond_trajectory_test_data.sh"
fi

echo "Test data for dimred/phate:"
echo "  input:  ${BEYOND_DIR}/proportions_output.h5mu  (obsm_input=proportions)"
echo ""
echo "Run the component manually:"
echo ""
echo "  viash run src/dimred/phate/config.vsh.yaml -- \\"
echo "    --input       ${BEYOND_DIR}/proportions_output.h5mu \\"
echo "    --obsm_input  proportions \\"
echo "    --output      ${BEYOND_DIR}/phate_output.h5mu"
