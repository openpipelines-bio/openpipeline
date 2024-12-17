import subprocess
from pathlib import Path
import mudata as md

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

print("> Running command with folder", flush=True)
input = meta["resources_dir"] + "/cellranger_tiny_fastq/htseq_counts.tsv"
reference = (
    meta["resources_dir"]
    + "/cellranger_tiny_fastq/cellranger_tiny_ref/genes/genes.gtf.gz"
)
output = "test_output.h5mu"

cmd_pars = [
    meta["executable"],
    "--input_id",
    "foo;bar",
    "--input_counts",
    f"{input};{input}",
    "--reference",
    reference,
    "--output",
    output,
    "---cpus",
    "2",
    "--output_compression",
    "gzip",
]
subprocess.run(cmd_pars, check=True)

print("> Check if file exists", flush=True)
output_path = Path(output)
assert output_path.is_file()

print("> Check contents", flush=True)
mdata = md.read_h5mu(output)

print(mdata)

assert "rna" in mdata.mod
assert mdata.n_obs == 2
assert mdata.mod["rna"].n_vars > 100

print("> Completed Successfully!", flush=True)
