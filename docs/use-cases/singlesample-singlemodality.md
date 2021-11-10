# Single sample - single modality

## RNA velocity

### Mapping with CellRanger

```
name: mapping/cellranger
container: litd/docker-cellranger:v6.1.1
input:
  fastqs: folder with fastq files
  transcriptome: folder with cellranger reference
output:
  output: folder with cellranger output
```

Example command
```bash
viash run src/1_cellranger/config.vsh.yaml -- \
  --input poc/pbmc_1k_v3_fastqs \
  --transcriptome poc/refdata-gex-GRCh38-2020-A \
  --output output/poc_output/ \
  --log testlog.txt
```

### Split into spliced/unspliced with velocyto

```
name: mapping/velocyto
container: python:3.8
python packages: [ velocyto ]
input:
  input: folder with cellranger output
  transcriptome: folder with cellranger reference
output:
  output: velocyto loom
```

Example command
```bash
viash run src/2_velocyto/config.vsh.yaml -- \
  --input output/poc_output/ \
  --transcriptome poc/refdata-gex-GRCh38-2020-A \
  --output output/poc_output.loom
```

### Compute RNA velocity with scvelo

```
name: mapping/velocyto
container: python:3.8
python packages: [ scvelo ]
input:
  input: velocyto loom
output:
  output: muon
    mod['velocity']:
      - layers['spliced']: Count matrix of spliced reads
      - layers['unspliced']: Count matrix of unspliced reads
      - layers['velocity']: Matrix of velocity vector
      ... more?
```

Example command
```bash
viash run src/3_scvelo/config.vsh.yaml -- \
  --input output/poc_output.loom \
  --output output/poc_output.h5mu
```