#!/bin/bash



# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_VER=21.10.6

# Add id to query .obs column with same name as .obs column of reference dataset containing batch labels
# Avoids nan values, which are not supported during harmony integration
viash run src/metadata/add_id/config.vsh.yaml -- \
--input resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu \
--output resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms_w_sample_id.h5mu \
--input_id pbmc_1k_protein_v3_mss \
--obs_output donor_assay

viash ns build -q harmony_knn

nextflow \
  run . \
  -main-script src/workflows/annotation/harmony_knn/test.nf \
  -entry test_wf \
  -resume \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config \
  -with-trace work/trace.txt