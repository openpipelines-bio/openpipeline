#!/bin/bash

## VIASH START
meta_executable="bin/viash run src/reference/make_reference/config.vsh.yaml --"
meta_resources_dir="resources_test"
## VIASH END

# create temporary directory
# tmpdir=$(mktemp -d "$meta_temp_dir/$meta_functionality_name-XXXXXXXX")
# function clean_up {
#     rm -rf "$tmpdir"
# }
# trap clean_up EXIT
tmpdir="."

start=30000
end=31500
seqkit subseq -r "$start:$end" "$meta_resources_dir/reference_gencodev41_chr1/reference.fa.gz" | gzip > "$tmpdir/reference_small.fa.gz"
zcat "$meta_resources_dir/reference_gencodev41_chr1/reference.gtf.gz" | \
  awk -v FS='\t' -v OFS='\t' "
    \$4 >= $start && \$5 <= $end {
      \$4 = \$4 - $start + 1;
      \$5 = \$5 - $start + 1;
      print;
    }" | \
  gzip > "$tmpdir/reference_small.gtf.gz"

echo "> Running $meta_functionality_name, writing to $tmpdir."
$meta_executable \
  --genome_fasta "$tmpdir/reference_small.fa.gz" \
  --gtf "$tmpdir/reference_small.gtf.gz" \
  --reference_archive "$tmpdir/myreference.tar.gz" \
  --extra_star_params "--genomeSAindexNbases 6" \
  ---cpus 2

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

echo ">> Checking whether output can be found"
[[ ! -f "$tmpdir/myreference.tar.gz" ]] && echo "Output tar file could not be found!" && exit 1

echo "> Test succeeded!"