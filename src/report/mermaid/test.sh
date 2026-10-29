#!/bin/bash

set -e

par_input="network.input"
par_output="output.png"

echo "Creating input file"
cat > "$par_input" << HERE
graph LR;
    A--> B & C & D;
    B--> A & E;
    C--> A & E;
    D--> A & E;
    E--> B & C & D;
HERE

echo "Running command"
"$meta_executable" --input "$par_input" --output "$par_output"

if [[ ! -f "$par_output" ]]; then
    echo "Output does not exist"
    exit 1
fi

echo "Checking that the output is a PNG file"
if ! head -c 8 "$par_output" | grep -q "PNG"; then
    echo "Output is not a PNG file"
    exit 1
fi

echo "Running command with --output_format svg (overriding the .png extension)"
"$meta_executable" --input "$par_input" --output "output_svg.png" --output_format "svg"

if [[ ! -f "output_svg.png" ]]; then
    echo "Output with --output_format svg does not exist"
    exit 1
fi

if ! grep -q "<svg" "output_svg.png"; then
    echo "Output with --output_format svg is not an SVG file"
    exit 1
fi

echo "Running command with --output_format pdf (overriding the .png extension)"
"$meta_executable" --input "$par_input" --output "output_pdf.png" --output_format "pdf"

if [[ ! -f "output_pdf.png" ]]; then
    echo "Output with --output_format pdf does not exist"
    exit 1
fi

if [[ "$(head -c 5 "output_pdf.png")" != "%PDF-" ]]; then
    echo "Output with --output_format pdf is not a PDF file"
    exit 1
fi

echo "Test succeeded!"