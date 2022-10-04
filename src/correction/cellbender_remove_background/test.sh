
#!/bin/bash

echo "> Running CellBender"

$meta_executable \
  --input "$meta_resources_dir/raw_feature_bc_matrix.h5" \
  --output "foo.h5" \
  --output_report "report.pdf" \
  --output_cell_barcodes "barcodes.csv" \
  --output_filtered "filtered.h5" \
  --epochs 5

[ -f foo.h5 ] || { echo "No foo.h5 output file found!"; exit 1; }
[ -f report.pdf ] || { echo "No report.pdf output file found!"; exit 1; }
[ -f barcodes.csv ] || { echo "No barcodes.csv output file found!"; exit 1; }
[ -f filtered.h5 ] || { echo "No filtered.h5 output file found!"; exit 1; }


echo "> Test succeeded!"