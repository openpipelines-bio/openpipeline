#!/bin/bash

nextflow run . \
  -main-script workflows/2_unimodal_singlesample/rna/main.nf \
  -profile docker \
  -resume \
  -entry test_wf