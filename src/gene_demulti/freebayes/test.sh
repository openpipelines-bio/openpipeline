set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --bam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --fasta_reference "$meta_resources_dir/demuxafy_test_data/refdata-cellranger-GRCh38-3.0.0/genome.fa" \
    --output freebayes_result/ \
    --vcf mixed.variant.vcf

[[ ! -f freebayes_result/mixed.variant.vcf ]] && echo "Output VCF file could not be found!" && exit 1

echo ">>> Test finished successfully"
