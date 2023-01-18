set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --bam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --bam_index "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam.bai" \
    --fasta "$meta_resources_dir/demuxafy_test_data/genome.fa" \
    --barcodes "$meta_resources_dir/demuxafy_test_data/barcodes.tsv" \
    --common_variants "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf" \
    --clusters 14  --threads 3 --skip_remap \
    --output soup_result/

[[ ! -f soup_result/clusters.tsv ]] && echo "Output donor assignment file could not be found!" && exit 1
[[ ! -f soup_result/assignment.tsv ]] && echo "Output summarized donor assignment file could not be found!" && exit 1

echo ">>> Test finished successfully"
