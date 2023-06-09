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
    exit
fi

echo "Test succeeded!"