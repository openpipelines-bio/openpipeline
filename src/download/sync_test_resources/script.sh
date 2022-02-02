#!/bin/bash

# extra_params=( )

# if [ "$par_verbose" != "true" ]; then 
#   extra_params+=( "--quiet" )
# fi

aws s3 sync s3://openpipelines-data "$par_output" --no-sign-request #"${extra_params[@]}"