#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

nextflow \
  run . \
  -main-script src/workflows/rna/rna_multisample/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/integration_tests.config

# Requires a CUDA-capable NVIDIA GPU. Listed under 'gpu_tests' in _viash.yaml, so the
# GitHub Actions integration test skips it; the Viash Hub CI has a GPU and runs it.
nextflow \
  run . \
  -main-script src/workflows/rna/rna_multisample/test.nf \
  -entry test_gpu_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/gpu.config \
  -c src/workflows/utils/integration_tests.config
