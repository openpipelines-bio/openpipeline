#!/bin/bash

## VIASH START
## VIASH END

cat > _viash2.yaml << EOM
info:
  test_resources:
    - type: s3
      path: s3://openpipelines-data/pbmc_1k_protein_v3
      dest: bar
EOM

echo ">> Run aws s3 sync"
"$meta_executable" \
  --input _viash2.yaml \
  --output foo \
  --exclude '*.h5' \
  --exclude '*.h5mu' \
  --exclude '*.mtx.gz' \
  --quiet

echo ">> Check whether the right files were copied"
[ ! -f foo/bar/pbmc_1k_protein_v3_metrics_summary.csv ] && echo csv should have been copied && exit 1
[ ! -f foo/bar/pbmc_1k_protein_v3_filtered_feature_bc_matrix/barcodes.tsv.gz ] && echo barcodes.tsv.gz should have been copied && exit 1
[ ! -f foo/bar/pbmc_1k_protein_v3_filtered_feature_bc_matrix/features.tsv.gz ] && echo features.tsv.gz should have been copied && exit 1
[ -f foo/bar/pbmc_1k_protein_v3_filtered_feature_bc_matrix/matrix.mtx.gz ] && echo matrix.mtx.gz should have been excluded && exit 1
[ -f foo/bar/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 ] && echo h5 should have been excluded && exit 1
[ -f foo/bar/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu ] && echo h5mu should have been excluded && exit 1

echo ">> Check whether content was found"
if ! grep -q 'Estimated Number of Cells' foo/bar/pbmc_1k_protein_v3_metrics_summary.csv; then
  echo Could not find content in csv
  exit 1
fi

echo ">> Test succeeded!"
