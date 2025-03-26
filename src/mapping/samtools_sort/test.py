import subprocess
from pathlib import Path

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

print("> sort and index", flush=True)
input = meta["resources_dir"] + "/cellranger_tiny_fastq/bam/possorted_genome_bam.bam"
output = "test_output.bam"

cmd_pars = [
    meta["executable"],
    "--input",
    input,
    "--output_bam",
    output,
    "--output_bai",
    output + ".bai",
    "---memory",
    "3gb",
    "---cpus",
    "2",
]
subprocess.run(cmd_pars, check=True)

print("> Check if file exists", flush=True)
bam_path = Path(output)
assert bam_path.is_file()
bai_path = Path(output + ".bai")
assert bai_path.is_file()

print("> sort by name (no index)", flush=True)
output2 = "test_output2.bam"

cmd_pars = [
    meta["executable"],
    "--input",
    input,
    "--output_bam",
    output2,
    "--sort_by_read_names",
    "---memory",
    "3gb",
    "---cpus",
    "2",
]
subprocess.run(cmd_pars, check=True)

print("> Check if file exists", flush=True)
bam_path2 = Path(output2)
assert bam_path2.is_file()

print("> Completed Successfully!", flush=True)
