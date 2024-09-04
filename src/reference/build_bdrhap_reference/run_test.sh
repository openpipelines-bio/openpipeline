#!/bin/bash



## VIASH START
meta_executable="bin/viash run src/reference/make_reference/config.vsh.yaml --"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

seqkit subseq -r 1:50000 "$meta_resources_dir/reference.fa.gz" | gzip > "$tmpdir/reference_small.fa.gz"
zcat "$meta_resources_dir/reference.gtf.gz" | awk '$4 < 50001 {print ;}' | gzip > "$tmpdir/reference_small.gtf.gz"


echo "> Running $meta_name, writing to $tmpdir."
$meta_executable \
  --genome_fasta "$tmpdir/reference_small.fa.gz" \
  --gtf "$tmpdir/reference_small.gtf.gz" \
  --reference_archive "$tmpdir/myreference.tar.gz" \
  ---cpus 2

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f "$tmpdir/myreference.tar.gz" ]] && echo "Output tar file could not be found!" && exit 1

echo "> Test succeeded!"