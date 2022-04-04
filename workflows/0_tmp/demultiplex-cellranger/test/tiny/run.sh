rm -rf out
mkdir out
mkdir -p work/bcl
tar -xvzf "$(pwd)/resources_test/cellranger-tiny-bcl-1.2.0/cellranger-tiny-bcl-1.2.0.tar.gz" -C "$(pwd)/work/bcl"

NXF_VER=20.12.1-edge ./bin/nextflow run ./workflows/demultiplex-cellranger/demultiplex-cellranger.nf \
  --input "$(pwd)/work/bcl/cellranger-tiny-bcl-1.2.0/" \
  --samplesheet "$(pwd)/resources_test/cellranger-tiny-bcl-1.2.0/cellranger-tiny-bcl-simple-1.2.0.csv" \
  --output "$(pwd)/out" \
  -resume \
  --mkfastq__memory "auto" \
  --mkfastq__cores "auto" -c local_nxf.config

./target/docker/test/integration_test/integration_test --test "./workflows/demultiplex-cellranger/test/tiny/test.py" --input "$(pwd)/out" --output "$(pwd)/out"
