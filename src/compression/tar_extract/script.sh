#!/usr/bin/env bash

## VIASH START

par_input='temp/test.tar'
par_output='temp/output_folder'
par_strip_components='0'
par_exclude='docs/figures'

## VIASH END

extra_params=()
mkdir -p $par_output # Create output directory if it doesn't exist already

if [ "$par_strip_components" != "" ]; then
    extra_params+=("--strip-components=$par_strip_components")
fi

if [ "$par_exclude" != "" ]; then
    extra_params+=("--exclude=$par_exclude")
fi

echo "Extracting $par_input to $par_output..."
echo ""
tar "${extra_params[@]}" -xvf "$par_input" -C "$par_output"
