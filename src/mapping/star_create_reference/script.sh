set -eo pipefail

extra_params=()

if [ ! -z "$sjdbGTFfile" ]; then
  extra_params+=("--sjdbGTFfile" "$sjdbGTFfile")
fi

if [ -z "$meta_cpus" ]; then
  meta_cpus=1
fi

STAR --runThreadN "$meta_cpus" \
    --runMode genomeGenerate \
    --genomeDir "$par_output" \
    --genomeFastaFiles "$par_input" \
    --genomeSAindexNbases "$par_genomeSAindexNbases" \
    "${extra_params[@]}"

