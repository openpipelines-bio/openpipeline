#!/bin/bash

nextflow run . \
  -main-script workflows/process_rna/multisample/main.nf \
  -profile docker \
  -resume \
  -entry test_wf