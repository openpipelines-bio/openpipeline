set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --input "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --vcf "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf" \
    --outputFile out.best.gt \
    --field GP 

[[ ! -f out.best.gt ]] && echo "Output result file could not be found!" && exit 1

echo ">>> Test finished successfully"