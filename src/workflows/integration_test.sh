#!/bin/bash

viash ns build --parallel --runner nextflow

Rscript src/workflows/integration_test.R | tee .viash_log_integration.txt