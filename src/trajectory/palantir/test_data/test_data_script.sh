#!/usr/bin/env bash
# Generate minimal test data for trajectory/palantir.
# Re-uses resources_test/beyond_test_data/atlas.h5mu (created by
# resources_test_scripts/beyond_trajectory_test_data.sh) — no new h5mu is needed.
#
# Usage:
#   cd /path/to/openpipeline
#   bash src/trajectory/palantir/test_data/test_data_script.sh
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
ATLAS="${REPO_ROOT}/resources_test/beyond_test_data/atlas.h5mu"

if [[ ! -f "${ATLAS}" ]]; then
  echo "ERROR: atlas.h5mu not found at ${ATLAS}"
  echo "Run resources_test_scripts/beyond_trajectory_test_data.sh first."
  exit 1
fi

echo "Test data for trajectory/palantir uses ${ATLAS}"
echo "No additional test files need to be generated."
