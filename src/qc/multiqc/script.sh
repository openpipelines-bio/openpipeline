#!/bin/bash

set -eo pipefail

multiqc -o "$par_output" $par_input
