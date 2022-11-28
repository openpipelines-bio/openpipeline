set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --bam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --output samtools_result/

[[ ! -f samtools_result/sorted.bam ]] && echo "Output preprocessed BAM file could not be found!" && exit 1

echo ">>> Test finished successfully"
