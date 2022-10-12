### VIASH START

par = {
par_input = "/path/to/input/bam_file.bam"
par_vcf = "/path/to/input/vcf_file.vcf"
par_numberOfSamples = 12 # The number of samples in the bam
par_sampleOutputFile = "/path/to/output_clustering_file.gz"
par_vcfOutputFile = "/path/to/output/cluster_associated_vcf_file.vcf.gz"
### VIASH END
set -x 

popscle dsc-pileup --sam "$par_input" --vcf "$par_vcf" --out pileup
popscle freemuxlet --plp pileup --out out --nsample "$par_numberOfSamples" 

cp out.clust1.samples.gz $par_sampleOutputFile
cp out.clust1.vcf.gz $par_vcfOutputFile
