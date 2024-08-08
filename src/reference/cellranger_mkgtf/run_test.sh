#!/bin/bash

set -eou pipefail

## VIASH START
meta_executable="bin/viash run src/reference/cellranger_mkgtf/config.vsh.yaml --"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

zcat "$meta_resources_dir/reference_gencodev41_chr1/reference.gtf.gz" | awk '$4 < 50001 {print ;}' | gzip > "$tmpdir/reference_small.gtf.gz"

expected_gene_types=("transcribed_unprocessed_pseudogene" "miRNA")
attribute_values=$(printf 'gene_type:%s,' "${expected_gene_types[@]}")
attribute_values=${attribute_values%,}  # remove trailing comma
echo $attribute_values

echo "> Running $meta_name, writing to $tmpdir."
$meta_executable \
  --input_gtf "$tmpdir/reference_small.gtf.gz" \
  --output_gtf "$tmpdir/myreference_filtered.gtf.gz" \
  --attribute "$attribute_values" \
  ---cpus ${meta_memory_gb:-1} \
  ---memory ${meta_memory_gb:-2}GB

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f "$tmpdir/myreference_filtered.gtf.gz" ]] && echo "Output gtf file could not be found!" && exit 1

echo ">> Checking attribute 'gene_type' in output gtf file"
unique_gene_types=$(zcat "$tmpdir/myreference_filtered.gtf.gz" | awk -F'\t' '$9 ~ /gene_type/ { split($9, a, ";"); for(i in a) if(a[i] ~ /gene_type/) print a[i] }' | sed 's/.*gene_type "\(.*\)".*/\1/' | sort -u)
echo "Expected gene types: ${expected_gene_types[@]}"
echo "Unique gene types: $unique_gene_types"
if [[ "${#expected_gene_types[@]}" != "$(echo "$unique_gene_types" | wc -w)" ]]; then
  echo "Error: Not all expected gene types were found in the output gtf file"
  exit 1
fi

echo "> Test succeeded!"