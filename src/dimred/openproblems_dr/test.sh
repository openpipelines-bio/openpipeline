#!/bin/bash

# viash ns build --parallel -q openproblems --setup cb

nextflow run . \
  -main-script target/nextflow/dimred/openproblems_dr/main.nf \
  -profile docker \
  --input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu \
  --publish_dir output/foo
