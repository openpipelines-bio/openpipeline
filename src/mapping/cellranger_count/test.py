import subprocess
from os import path
import logging
from sys import stdout

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

logger.info("> Running command with folder")
input = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_fastq/"
reference = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_ref/"
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
assert path.exists(output + "/filtered_feature_bc_matrix.h5"), "No output was created."


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
assert path.exists(output + "/filtered_feature_bc_matrix.h5"), "No output was created."



logger.info("> Completed Successfully!")