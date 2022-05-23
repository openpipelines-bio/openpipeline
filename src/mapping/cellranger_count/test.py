import subprocess
from os import path

## VIASH START
meta = {
    "functionality_name": "cellranger_count",
    "resources_dir": "resources_test"
}
## VIASH END

print("> Running command")
input = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_fastq/"
reference = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_ref/"
output = "test_output"

cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--reference", reference,
    "--output", output,
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

print("> Check if file exists")
assert path.exists(output + "/filtered_feature_bc_matrix.h5"), "No output was created."

print("> Completed Successfully!")