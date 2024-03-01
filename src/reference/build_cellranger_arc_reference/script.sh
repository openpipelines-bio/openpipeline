#!/bin/bash

set -eo pipefail

## VIASH START
par_genome_fasta="resources_test/reference_gencodev41_chr1/reference.fa.gz"
par_annotation_gtf="resources_test/reference_gencodev41_chr1/reference.gtf.gz"
par_output="gencode_v41_annotation_cellranger.tar.gz"
par_motifs_file="resources_test/reference_gencodev41_chr1/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt.modified"
## VIASH END

# create temporary directory
tmpdir=$(mktemp -d "$VIASH_TEMP/$meta_functionality_name-XXXXXXXX")
function clean_up {
    rm -rf "$tmpdir"
}
trap clean_up EXIT

# just to make sure
echo "> Getting path of fasta file"
par_genome_fasta=`realpath $par_genome_fasta`
echo "> Getting path of annotation file"
par_annotation_gtf=`realpath $par_annotation_gtf`
echo "> Getting path of output file"
par_output=`realpath $par_output`
echo "> Getting path of motifs file"
par_motifs_file=`realpath $par_motifs_file`

# process params
extra_params=( )

if [ ! -z "$meta_cpus" ]; then 
  extra_params+=( "nthreads: \"$meta_cpus"\" )
fi
if [ ! -z "$meta_memory_gb" ]; then 
  # always keep 2gb for the OS itself
  memory_gb=`python -c "print(int('$meta_memory_gb') - 2)"`
  extra_params+=( "memgb: \"$memory_gb"\" )
fi

echo "> Unzipping input files"
unpigz -c "$par_genome_fasta" > "$tmpdir/genome.fa"

echo "> Building star index"
cd "$tmpdir"

echo "> Building config"
config_in="${tmpdir}/config"
echo """{
    organism: \"Homo_sapiens\"
    genome: [\"GRCh38\"]
    input_fasta: [\""${tmpdir}/genome.fa"\"]
    input_gtf: [\""${par_annotation_gtf}\""]
    non_nuclear_contigs: [\""chrM\""]
    input_motifs: \""$par_motifs_file"\"
    $(printf "%s\n" "${extra_params[@]}")
}""" > "$config_in"

echo "> Config content:"
cat ${config_in}

echo "> Running cellranger"
cellranger-arc mkref --config=${config_in}

echo "> Creating archive"
tar --use-compress-program="pigz -k " -cf "$par_output" -C "${tmpdir}/GRCh38" .

