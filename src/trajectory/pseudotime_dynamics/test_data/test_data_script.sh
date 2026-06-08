#!/usr/bin/env bash
# Test data for trajectory/pseudotime_dynamics.
# Delegates to resources_test_scripts/beyond_trajectory_test_data.sh.
# Input: resources_test/beyond_test_data/pseudotime_output.h5mu
#   (has obs["palantir_pseudotime"] and uns["proportions"],
#    as produced by trajectory/palantir after metadata/calculate_proportions)
#
# Run from the repo root:
#   bash src/trajectory/pseudotime_dynamics/test_data/test_data_script.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
BEYOND_DIR="${REPO_ROOT}/resources_test/beyond_test_data"

if [[ ! -f "${BEYOND_DIR}/pseudotime_output.h5mu" ]]; then
  echo "Generating BEYOND test data..."
  bash "${REPO_ROOT}/resources_test_scripts/beyond_trajectory_test_data.sh"
fi

echo "Test data for trajectory/pseudotime_dynamics:"
echo "  input:  ${BEYOND_DIR}/pseudotime_output.h5mu"
echo ""
echo "Run the component manually:"
echo ""
echo "  viash run src/trajectory/pseudotime_dynamics/config.vsh.yaml -- \\"
echo "    --input  ${BEYOND_DIR}/pseudotime_output.h5mu \\"
echo "    --output ${BEYOND_DIR}/dynamics_output.h5mu"
