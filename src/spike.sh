#!/bin/bash

viash run src/fetch/download_file/config.vsh.yaml -- \
  --input https://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_protein_v3_nextgem/5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5 \
  --output resources/muon_spike/5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5

viash run src/convert/h510x_to_h5mu/config.vsh.yaml -- \
  --input resources/muon_spike/5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5 \
  --output resources/muon_spike/5k_pbmc_protein_v3_nextgem_filtered_feature_bc_matrix.h5mu