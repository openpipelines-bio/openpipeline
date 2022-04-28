#!/bin/bash

## VIASH START
par_input='s3://openpipelines-data'
par_output='resources_test'
## VIASH END

# extra_params=( )

# if [ "$par_verbose" != "true" ]; then
#   extra_params+=( "--quiet" )
# fi

aws s3 sync "$par_input" "$par_output" --no-sign-request #"${extra_params[@]}"
