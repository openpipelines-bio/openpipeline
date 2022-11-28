set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --samFile "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --barcodeFile "$meta_resources_dir/demuxafy_test_data/barcodes.tsv" \
    --output cellsnp_result/ \
    --regionsVCF "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf"

[[ ! -f cellsnp_result/cellSNP.base.vcf ]] && echo "Output VCF file could not be found!" && exit 1

echo ">>> Test finished successfully"
