#!/usr/bin/env bash
# Test data for cluster/cellular_communities.
# Delegates to resources_test_scripts/beyond_test_data.sh.
# Input: resources_test/beyond_test_data/dynamics_output.h5mu
#   (has uns["proportions"] and uns["dynamics"],
#    as produced by trajectory/pseudotime_dynamics)
#
# Run from the repo root:
#   bash src/cluster/cellular_communities/test_data/test_data_script.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
BEYOND_DIR="${REPO_ROOT}/resources_test/beyond_test_data"

if [[ ! -f "${BEYOND_DIR}/dynamics_output.h5mu" ]]; then
  echo "Generating BEYOND test data..."
  bash "${REPO_ROOT}/resources_test_scripts/beyond_test_data.sh"
fi

echo "Test data for cluster/cellular_communities:"
echo "  input:  ${BEYOND_DIR}/dynamics_output.h5mu"
echo ""
echo "Run the component manually:"
echo ""
echo "  viash run src/cluster/cellular_communities/config.vsh.yaml -- \\"
echo "    --input         ${BEYOND_DIR}/dynamics_output.h5mu \\"
echo "    --n_communities 3 \\"
echo "    --output        ${BEYOND_DIR}/communities_output.h5mu"
