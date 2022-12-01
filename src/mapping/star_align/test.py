import subprocess
from os import path
import logging
from sys import stdout
from pathlib import Path

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)


## VIASH START
meta = {
    "functionality_name": "cellranger_count",
    "resources_dir": "resources_test"
}
## VIASH END

## Test 1: use input dir
logger.info("> Running command with folder")
input = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_fastq/"
reference = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_ref/star/"
output = "test_output"

cmd_pars = [
    meta["executable"],
    "--input", input,
    "--reference", reference,
    "--output", output,
    "---cores", "2",
    "---memory", "5gb"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")

output_path = Path(output)
assert (output_path / "Log.final.out" ).is_file(), "No output log was created."
assert (output_path / "SJ.out.tab" ).is_file(), "No output file was created."


## Test 2: use input files
logger.info("> Running command with fastq files")
input_files = [
  input + "tinygex_S1_L001_R1_001.fastq.gz",
  input + "tinygex_S1_L001_R2_001.fastq.gz"
]
output = "test_output2"

cmd_pars = [
    meta["executable"],
    "--input", input_files[0],
    "--input", input_files[1],
    "--reference", reference,
    "--output", output,
    "---cores", "2",
    "---memory", "5gb"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")
assert path.exists(output + "/Log.final.out"), "No output log was created."
assert path.exists(output + "/SJ.out.tab"), "No output was created."


logger.info("> Completed Successfully!")