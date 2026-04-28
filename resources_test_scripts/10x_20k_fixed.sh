#!/bin/bash

set -eo pipefail

# ensure that the command below is run from the root of the repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# settings
ID=10x_5k_fixed
OUT="resources_test/$ID"

# create raw directory
raw_dir="$OUT/raw"
mkdir -p "$raw_dir"

# Check whether seqkit is available
if ! command -v seqkit &> /dev/null; then
    echo "This script requires seqkit. Please make sure the binary is added to your PATH."
    exit 1
fi

# check whether reference is available
reference_dir="resources_test/reference_gencodev41_chr1/"
genome_tar="$reference_dir/reference_cellranger.tar.gz"
if [[ ! -f "$genome_tar" ]]; then
    echo "$genome_tar does not exist. Please create the reference genome first"
    exit 1
fi

# create tempdir
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}

# dataset page:
# https://www.10xgenomics.com/datasets/mixture-of-healthy-and-cancer-ffpe-tissues-dissociated-using-miltenyi-ffpe-tissue-dissociation-kit-multiplexed-samples-4-probe-barcodes-1-standard

# download and untar source fastq files
tar_dir="$HOME/.cache/openpipeline/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex"
if [[ ! -d "$tar_dir" ]]; then
    mkdir -p "$tar_dir"

    # download fastqs and untar
    wget "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/7.1.0/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_Multiplex/4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex_Multiplex_fastqs.tar" -O "$tar_dir.tar"
    tar -xvf "$tar_dir.tar" -C "$tar_dir" --strip-components=1
    rm "$tar_dir.tar"
fi

function seqkit_head {
  input="$1"
  output="$2"
  if [[ ! -f "$output" ]]; then
    echo "> Processing `basename $input`"
    seqkit head -n 200000 "$input" | gzip > "$output"
  fi
}

orig_sample_id="4plex_human_liver_colorectal_ovarian_panc_scFFPE_multiplex"

seqkit_head "$tar_dir/${orig_sample_id}_S1_L001_R1_001.fastq.gz" "$raw_dir/${orig_sample_id}_subset_S1_L001_R1_001.fastq.gz"
seqkit_head "$tar_dir/${orig_sample_id}_S1_L001_R2_001.fastq.gz" "$raw_dir/${orig_sample_id}_subset_S1_L001_R2_001.fastq.gz"

# download feature reference
feature_ref="$raw_dir/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_feature_reference.csv"
if [[ ! -f "$feature_ref" ]]; then
  wget "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Multiplex/4plex_mouse_LymphNode_Spleen_TotalSeqC_multiplex_Multiplex_count_feature_reference.csv" -O "$feature_ref"
fi

# download probe set
probe_set="$raw_dir/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv"
if [[ ! -f "$probe_set" ]]; then
  wget "https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv" -O "$probe_set"
fi

sed -i 's/#reference_genome=GRCh38/#reference_genome=output/g' "$probe_set"

probe_set_corrected="$raw_dir/Chromium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A_corrected.csv"
if [[ ! -f "$probe_set_corrected" ]]; then
  reference_gtf="resources_test/reference_gencodev41_chr1/reference.gtf.gz"
  gunzip -c "$reference_gtf" > "$TMPDIR/uncompressed_ref.gtf" 
  cat "$probe_set" | while read line || [[ -n $line ]];
  do
    echo "Line: $line"
    old_id=$( printf "%s\n" "$line" | awk -F',' '{print $1}' )
    echo "Old ID: $old_id"
    if [[ "$old_id" == "gene_id" ]] || [[ "$old_id" == \#* ]] ; then
      echo "Just writing line"
      printf "%s\n" "$line" >> "$probe_set_corrected"
    else
      gtf_lookup=$(grep "$old_id" "$TMPDIR/uncompressed_ref.gtf" || test $? = 1;)
      if [ ! -z "$gtf_lookup" ]; then
        echo "Found hit"
        new_id=$(echo "$gtf_lookup" | awk '{if ($3 == "gene") print $10;}' | sed -e "s/^\"//" -e "s/\";$//")
        echo "New ID: $new_id"
        new_line=${line/"$old_id"/"$new_id"}
        printf "%s\n" "$new_line" >> "$probe_set_corrected"
      else
        echo "Did not find hit"
      fi
    fi
  done
fi

# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
# Output FASTA:
#   >chr1 1
input_fastq="$HOME/.cache/openpipeline/GRCh38.primary_assembly.genome.fa.gz"
fasta_modified="$TMPDIR/GRCh38.primary_assembly.genome.modified.fa"
if [[ ! -f "$input_fastq" ]]; then
  wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/GRCh38.primary_assembly.genome.fa.gz" -O "$input_fastq"
fi
zcat "$input_fastq" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

pigz --fast "$fasta_modified"
fasta_modified="$fasta_modified.gz"
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
input_gtf="$HOME/.cache/openpipeline/gencode.v41.annotation.gtf.gz"
gtf_modified="$TMPDIR/gencode.v41.annotation.gtf.modified.gtf"
if [[ ! -f "$input_gtf" ]]; then
  wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz" -O "$input_gtf"
fi

REGEX="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
zcat "$input_gtf" \
    | sed -E 's/gene_id "'"$REGEX"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$REGEX"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$REGEX"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"
pigz --fast "$gtf_modified"
gtf_modified="$gtf_modified.gz"

final_genome="$HOME/.cache/openpipeline/GRCh38.cellranger.genome.fa.gz"
if [ ! -f "$final_genome" ]; then
  nextflow \
    run openpipelines-bio/openpipeline \
    -r 3.0.0 \
    -main-script target/nextflow/workflows/ingestion/make_reference/main.nf \
    -profile docker \
    -resume \
    --id "GRCh38" \
    --genome_fasta "$fasta_modified" \
    --transcriptome_gtf "$gtf_modified" \
    --target "cellranger" \
    --output_fasta "reference.fa.gz" \
    --output_gtf "reference.gtf.gz" \
    --output_cellranger "GRCh38.cellranger.genome.fa.gz" \
    --publish_dir "$HOME/.cache/openpipeline/"
fi


# Run mapping pipeline
cat > /tmp/params.yaml << HERE
param_list:
- id: "$ID"
  input: "$raw_dir"
  library_id:
    - ${orig_sample_id}_subset
  library_type:
    - "Gene Expression"
  library_lanes:
    - "any"

probe_set: "$probe_set_corrected"
gex_reference: "$genome_tar"
feature_reference: "$feature_ref"
probe_barcode_ids:
  - BC001
  - BC002
  - BC003
  - BC004
sample_ids:
  - Liver_BC1
  - Ovarian_BC2
  - Colorectal_BC3
  - Pancreas_BC4
gex_generate_bam: false
sample_force_cells:
  - 5000
  - -1
  - -1
  - -1
HERE

nextflow \
  run openpipelines-bio/openpipeline \
  -r 3.0.0 \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "$OUT/processed"

mkdir -p "${OUT}/processed"

nextflow \
  run . \
  -main-script target/nextflow/mapping/cellranger_multi/main.nf \
  -resume \
  -profile docker,mount_temp \
  -params-file /tmp/params.yaml \
  -c src/workflows/utils/labels_ci.config \
  --publish_dir "${OUT}_v10/processed"
