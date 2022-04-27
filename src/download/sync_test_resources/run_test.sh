

./$meta_functionality_name \
  --input s3://openpipelines-data/pbmc_1k_protein_v3 \
  --output foo

[ ! -f foo/pbmc_1k_protein_v3_metrics_summary.csv ] && echo Output file could not be found && exit 1

if ! grep -q 'Estimated Number of Cells' foo/pbmc_1k_protein_v3_metrics_summary.csv; then
  echo Could not find content
  exit 1
fi

echo Test succeeded!