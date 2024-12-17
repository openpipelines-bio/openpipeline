import subprocess
from pathlib import Path
import pandas as pd


## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

print("> Running command with folder", flush=True)
input = meta["resources_dir"] + "/cellranger_tiny_fastq/bam/possorted_genome_bam.bam"
reference = (
    meta["resources_dir"]
    + "/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz"
)
output = "test_output.tsv"

cmd_pars = [
    meta["executable"],
    "--input",
    input,
    "--reference",
    reference,
    "--output",
    output,
    "---cpus",
    "2",
]
subprocess.run(cmd_pars, check=True)

print("> Check if file exists", flush=True)
output_path = Path(output)
assert output_path.is_file()

print("> Check contents", flush=True)
counts = pd.read_table(output_path, sep="\t")
assert counts.shape[0] > 100
assert counts.shape[1] == 2

print("> Completed Successfully!", flush=True)
