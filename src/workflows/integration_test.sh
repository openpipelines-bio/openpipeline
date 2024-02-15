#!/bin/bash

viash ns build --parallel -p nextflow

Rscript src/workflows/integration_test.R | tee .viash_log_integration.txt