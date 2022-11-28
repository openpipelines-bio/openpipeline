set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --bam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --vcf "$meta_resources_dir/demuxafy_test_data/mixed.variant.vcf" \
    --num 12 \
    --bar "$meta_resources_dir/demuxafy_test_data/barcodes.tsv" \
    --com "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf" \
    --ref ref_filtered.csv \
    --alt alt_filtered.csv \
    --output scSplit_result/

[[ ! -f scSplit_result/scSplit_result.csv ]] && echo "Output donor assignment file could not be found!" && exit 1

echo ">>> Test finished successfully"
