#!/bin/bash

# Test data for report/deseq2_report: output of differential_expression/deseq2
# run with --export_normalized_counts on the pseudobulk test fixture
# (resources_test/annotation_test_data/TS_Blood_filtered_pseudobulk.h5mu),
# restricted to the 500 most highly expressed genes that are detected in every
# sample, to keep the files small enough to keep in the repository.
#   overall/        one analysis over all samples (~ cell_type + treatment)
#   per_cell_type/  the erythrocyte results of an analysis per cell type
#                   (--obs_cell_group cell_type)
# The treatment labels of the fixture are random, so few or no genes are
# differentially expressed.

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

IN=resources_test/annotation_test_data/TS_Blood_filtered_pseudobulk.h5mu
OUT=src/report/deseq2_report/test_data
N_GENES=500

[ -f "$IN" ] || { echo "Missing $IN, sync the test resources first"; exit 1; }

TMP_DIR=$(mktemp -d "$PWD/tmp_deseq2_report_XXXXXX")
function clean_up {
  rm -rf "$TMP_DIR"
}
trap clean_up EXIT

echo "> Keeping the $N_GENES most highly expressed genes detected in every sample"
cp "$IN" "$TMP_DIR/pseudobulk.h5mu"
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$TMP_DIR:/data" -w /data python:3.12-slim bash -c "
pip install --quiet --user mudata > /dev/null 2>&1
python - <<'HEREDOC'
import mudata as mu
import numpy as np

mdata = mu.read_h5mu('pseudobulk.h5mu')
counts = np.asarray(mdata['rna'].X)
detected = (counts > 0).all(axis=0)
mean = np.where(detected, counts.mean(axis=0), -1)
keep = np.sort(np.argsort(mean)[::-1][:$N_GENES])
subset = mu.MuData({'rna': mdata['rna'][:, keep].copy()})
subset.write_h5mu('pseudobulk_subset.h5mu')
print(subset)
HEREDOC
"

echo "> Running DESeq2 on all samples"
rm -rf "$OUT/overall" "$OUT/per_cell_type"
viash run src/differential_expression/deseq2/config.vsh.yaml --engine docker -- \
  --input "$TMP_DIR/pseudobulk_subset.h5mu" \
  --design_formula "~ cell_type + treatment" \
  --contrast_column treatment \
  --contrast_values stim \
  --contrast_values ctrl \
  --output_dir "$OUT/overall" \
  --export_normalized_counts

echo "> Running DESeq2 per cell type, keeping the erythrocyte results"
viash run src/differential_expression/deseq2/config.vsh.yaml --engine docker -- \
  --input "$TMP_DIR/pseudobulk_subset.h5mu" \
  --obs_cell_group cell_type \
  --design_formula "~ treatment" \
  --contrast_column treatment \
  --contrast_values stim \
  --contrast_values ctrl \
  --output_dir "$TMP_DIR/per_cell_type" \
  --export_normalized_counts
mkdir -p "$OUT/per_cell_type"
cp "$TMP_DIR"/per_cell_type/deseq2_analysis_erythrocyte* "$OUT/per_cell_type/"

echo "> Done: $OUT"
