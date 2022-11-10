set -ex
echo "$meta_resources_dir"
echo ">>> Running executable"
$meta_executable \
    --input "$meta_resources_dir/demuxafy_test_data/pooled.sorted.bam" \
    --vcf "$meta_resources_dir/demuxafy_test_data/test_dataset.vcf" \
    --sampleOutputFile outputSamples.gz \
    --vcfOutputFile output.vcf.gz \
    --numberOfSamples 2 

[[ ! -f outputSamples.gz ]] && echo "Output sample file could not be found!" && exit 1
[[ ! -f output.vcf.gz ]] && echo "Output VCF file could not be found!" && exit 1

echo ">>> Test finished successfully"