#!/bin/bash



## VIASH START
meta_executable="bin/viash run src/reference/make_reference/config.vsh.yaml --"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

function seqkit_head {
  input="$1"
  output="$2"
  if [[ ! -f "$output" ]]; then
    echo "> Processing `basename $input`"
    seqkit subseq -r 1:50000 "$input" | gzip > "$output"
  fi
}

seqkit_head "$meta_resources_dir/reference_gencodev41_chr1/reference.fa.gz" "$tmpdir/reference_small.fa.gz"
zcat "$meta_resources_dir/reference_gencodev41_chr1/reference.gtf.gz" | awk '$4 < 50001 {print ;}' | gzip > "$tmpdir/reference_small.gtf.gz"

echo "> Running $meta_functionality_name, writing to $tmpdir."
$meta_executable \
  --genome_fasta "$tmpdir/reference_small.fa.gz" \
  --transcriptome_gtf "$tmpdir/reference_small.gtf.gz" \
  --output "$tmpdir/myreference.tar.gz" \
  ---cpus ${meta_memory_gb:-1} \
  ---memory ${meta_memory_gb:-5}GB

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f "$tmpdir/myreference.tar.gz" ]] && echo "Output tar file could not be found!" && exit 1

echo "> Test succeeded!"