set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --sam "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --vcf "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf" \
    --output freemuxlet_result/ \
    --outDsc freemux_out_samples \
    --out freemx_out \
    --nsample 2

[[ ! -f freemuxlet_result/freemux_out_samples.plp.gz ]] && echo "Output sample file could not be found!" && exit 1
[[ ! -f freemuxlet_result/freemux_out.clust1.samples.gz ]] && echo "Output VCF file could not be found!" && exit 1

echo ">>> Test finished successfully"
