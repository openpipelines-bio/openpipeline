#!/usr/bin/env bash
# Test data for stats/trait_associations.
# Delegates to resources_test_scripts/beyond_trajectory_test_data.sh.
# Inputs:
#   resources_test/beyond_test_data/proportions_output.h5mu — has uns["proportions"]
#   resources_test/beyond_test_data/traits.csv              — 8 donors × 5 traits
#
# Run from the repo root:
#   bash src/stats/trait_associations/test_data/test_data_script.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
BEYOND_DIR="${REPO_ROOT}/resources_test/beyond_test_data"

if [[ ! -f "${BEYOND_DIR}/proportions_output.h5mu" ]]; then
  echo "Generating BEYOND test data..."
  bash "${REPO_ROOT}/resources_test_scripts/beyond_trajectory_test_data.sh"
fi

echo "Test data for stats/trait_associations:"
echo "  input:      ${BEYOND_DIR}/proportions_output.h5mu"
echo "  traits_csv: ${BEYOND_DIR}/traits.csv"
echo ""
echo "Run the component manually:"
echo ""
echo "  viash run src/stats/trait_associations/config.vsh.yaml -- \\"
echo "    --input         ${BEYOND_DIR}/proportions_output.h5mu \\"
echo "    --traits_csv    ${BEYOND_DIR}/traits.csv \\"
echo "    --trait_columns age --trait_columns pmi \\"
echo "    --output        ${BEYOND_DIR}/trait_associations_output.h5mu \\"
echo "    --output_csv    ${BEYOND_DIR}/trait_associations.csv"
