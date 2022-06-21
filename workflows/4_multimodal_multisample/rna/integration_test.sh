#!/bin/bash

nextflow run . \
  -main-script workflows/4_multimodal_multisample/rna/main.nf \
  -profile docker \
  -resume \
  -entry test_wf