import subprocess
from os import path
import sys

## VIASH START
meta = {
    "functionality_name": "cellranger_count",
    "resources_dir": "resources_test"
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

logger.info("> Running command with folder")
input = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_fastq/"
reference = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_ref/"
output = "test_output"

cmd_pars = [
    meta["executable"],
    "--input", input,
    "--reference", reference,
    "--output", output
]
if meta.get("cpus"):
    cmd_pars.extend(["---cpus", str(meta["cpus"])])
if meta.get("memory_gb"):
    cmd_pars.extend(["---memory", f"{meta['memory_gb']}gb"])

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
    "--output", output
]
if meta.get("cpus"):
    cmd_pars.extend(["---cpus", str(meta["cpus"])])
if meta.get("memory_gb"):
    cmd_pars.extend(["---memory", f"{meta['memory_gb']}gb"])
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")
assert path.exists(output + "/filtered_feature_bc_matrix.h5"), "No output was created."



logger.info("> Completed Successfully!")