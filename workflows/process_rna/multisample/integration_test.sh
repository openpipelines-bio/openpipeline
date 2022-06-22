#!/bin/bash

nextflow run . \
  -main-script workflows/3_unimodal_multisample/rna/main.nf \
  -profile docker \
  -resume \
  -entry test_wf