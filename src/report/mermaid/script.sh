#!/bin/bash

mmdc -p "$meta_resources_dir/puppeteer-config.json" \
    -i "$par_input" \
    -o "$par_output" \
    --width "$par_width" \
    --height "$par_height" \
    ${par_background_color:+--backgroundColor $par_background_color} \
    ${output_format:+--outputFormat $par_output_format}

