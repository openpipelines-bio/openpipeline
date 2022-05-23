#!/bin/bash

aws s3 sync --profile di "resources_test" "s3://openpipelines-data" --delete --dryrun