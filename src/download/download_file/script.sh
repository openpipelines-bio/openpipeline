#!/bin/bash

extra_params=( )

if [ "$par_verbose" != "true" ]; then 
  extra_params+=( "--quiet" )
fi

wget "$par_input" -O "$par_output" "${extra_params[@]}"