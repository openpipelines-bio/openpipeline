#!/bin/bash

mmdc -p "$meta_resources_dir/puppeteer-config.json" \
    -i "$par_input" \
    -o "$par_output" \
    --size "$par_size" \
    ${par_background_color:+--backgroundColor $par_background_color} \
    ${par_output_format:+--outputFormat "$par_output_format"}

