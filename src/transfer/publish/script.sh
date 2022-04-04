#!/bin/bash

## VIASH START
par_input="input.txt"
par_output="output.txt"
## VIASH END

# We can have dirs in the output parameter
# Let's cope with those...

if [[ ! $(dirname "ddd/aaa") == "." ]]; then
  dirs=$(dirname $par_output)
  f=$(basename $par_output)
  mkdir -p $dirs
else
  f=$par_output
fi

cp "$par_input" "$par_output"
