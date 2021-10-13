#!/bin/bash

../bin/viash build ../src/report/mermaid/config.vsh.yaml -o ../target/report/mermaid/

rm -rf figures
mkdir figures

for figure in $(ls -1 ./mermaid/ | grep .mermaid | sed 's/\.[^.]*$//'); do
    echo "Generating mermaid figure ${figure}"
    ../target/report/mermaid/mermaid -i ./mermaid/${figure}.mermaid -o ./figures/${figure}.png
done
