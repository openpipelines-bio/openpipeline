#!/bin/bash

## VIASH START
par_input='s3://openpipelines-data'
par_output='resources_test'
## VIASH END

# extra_params=( )

# if [ "$par_verbose" != "true" ]; then
#   extra_params+=( "--quiet" )
# fi

# Disable the use of the Amazon EC2 instance metadata service (IMDS).
# see https://florian.ec/blog/github-actions-awscli-errors/
# or https://github.com/aws/aws-cli/issues/5234#issuecomment-705831465
export AWS_EC2_METADATA_DISABLED=true

aws s3 sync "$par_input" "$par_output" --no-sign-request #"${extra_params[@]}"
