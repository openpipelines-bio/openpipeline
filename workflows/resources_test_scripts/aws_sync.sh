#!/bin/bash

aws s3 sync --profile di "resources_test" "s3://openpipelines-data" --exclude */temp_* --delete --dryrun