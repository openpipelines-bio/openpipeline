#!/bin/bash

set -eo pipefail

aws s3 sync --profile di "resources_test" "s3://openpipelines-data" --exclude "temp_*" --exclude "tmp_*" --delete --dryrun

id=cellranger_tiny_fastq
aws s3 sync --profile di "resources_test/$id" "s3://openpipelines-data/$id" --exclude "temp_*" --exclude "tmp_*" --delete --dryrun