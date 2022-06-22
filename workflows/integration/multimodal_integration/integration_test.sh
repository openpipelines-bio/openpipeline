#!/bin/bash

nextflow run . \
  -main-script workflows/integration/multimodal_integration/main.nf \
  -profile docker \
  -resume \
  -entry test_wf