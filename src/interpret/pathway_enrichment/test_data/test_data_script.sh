#!/usr/bin/env bash
# Test data for interpret/pathway_enrichment.
# Delegates to resources_test_scripts/beyond_trajectory_test_data.sh.
# Inputs:
#   resources_test/beyond_test_data/atlas.h5mu      — h5mu with gene names GENE00000–GENE01999
#   resources_test/beyond_test_data/de_ExN.csv      — DESeq2 results (same gene names)
#   resources_test/beyond_test_data/gene_sets.gmt   — 6 gene sets (SIG_A/B/C/ExN/InN/Ast)
#
# Run from the repo root:
#   bash src/interpret/pathway_enrichment/test_data/test_data_script.sh

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
BEYOND_DIR="${REPO_ROOT}/resources_test/beyond_test_data"

if [[ ! -f "${BEYOND_DIR}/atlas.h5mu" ]]; then
  echo "Generating BEYOND test data..."
  bash "${REPO_ROOT}/resources_test_scripts/beyond_trajectory_test_data.sh"
fi

echo "Test data for interpret/pathway_enrichment:"
echo "  h5mu:      ${BEYOND_DIR}/atlas.h5mu"
echo "  de genes:  ${BEYOND_DIR}/de_ExN.csv"
echo "  gene sets: ${BEYOND_DIR}/gene_sets.gmt"
echo ""
echo "Run the component manually:"
echo ""
echo "# prerank GSEA:"
echo "  viash run src/interpret/pathway_enrichment/config.vsh.yaml -- \\"
echo "    --input           ${BEYOND_DIR}/atlas.h5mu \\"
echo "    --input_degenes   ${BEYOND_DIR}/de_ExN.csv \\"
echo "    --method          prerank \\"
echo "    --gene_sets       ${BEYOND_DIR}/gene_sets.gmt \\"
echo "    --permutation_num 100 \\"
echo "    --output          ${BEYOND_DIR}/pathway_prerank_output.h5mu \\"
echo "    --output_csv_dir  ${BEYOND_DIR}/pathway_prerank_results/"
echo ""
echo "# ORA:"
echo "  viash run src/interpret/pathway_enrichment/config.vsh.yaml -- \\"
echo "    --input           ${BEYOND_DIR}/atlas.h5mu \\"
echo "    --input_degenes   ${BEYOND_DIR}/de_ExN.csv \\"
echo "    --method          ora \\"
echo "    --gene_sets       ${BEYOND_DIR}/gene_sets.gmt \\"
echo "    --output          ${BEYOND_DIR}/pathway_ora_output.h5mu \\"
echo "    --output_csv_dir  ${BEYOND_DIR}/pathway_ora_results/"
