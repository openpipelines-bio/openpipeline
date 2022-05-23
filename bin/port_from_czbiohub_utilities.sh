#!/bin/bash

tmpdir=$(mktemp -d)
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

curl -L -s "https://github.com/czbiohub/utilities/archive/refs/heads/main.zip" -o "$tmpdir/main.zip"
unzip -q "$tmpdir/main.zip" -d "$tmpdir"

dirs=(
    "src/mapping/cellranger_count/"
    "src/mapping/cellranger_count_split/"
    "src/demux/cellranger_mkfastq/"
    "workflows/1_ingestion/cellranger/"
    "workflows/1_ingestion/cellranger_demux/"
    "workflows/1_ingestion/cellranger_mapping/"
    "workflows/resources_test_scripts/cellranger_tiny_bcl.sh"
    "workflows/resources_test_scripts/cellranger_tiny_fastq.sh"
)

for dir in ${dirs[@]}; do
    echo "############## Syncing $dir ##############"
    mkdir -p `dirname "$dir"`
    rsync -avz --delete "$tmpdir/utilities-main/$dir" "$dir"
    echo
done