set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --bam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --fasta_reference "$meta_resources_dir/demuxafy_test_data/genome.fa" \
    --region 1 --output freebayes_result/ \
    --vcf mixed_variant_chr1.vcf

[[ ! -f freebayes_result/mixed_variant_chr1.vcf ]] && echo "Output VCF file could not be found!" && exit 1

echo ">>> Test finished successfully"
