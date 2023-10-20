import subprocess
from os import path
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
import shutil

## VIASH START
meta = {
    "functionality_name": "cellranger_count",
    "resources_dir": "resources_test"
}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

## Test 1: use input dir
logger.info("> Running command with folder")
input = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_fastq/"
input_files = [
  input + "tinygex_S1_L002_R1_001.fastq.gz",
  input + "tinygex_S1_L002_R2_001.fastq.gz"
]
reference = meta["resources_dir"] + "/cellranger_tiny_fastq/cellranger_tiny_ref_v2_7_10_a/"
output = "test_output"

with TemporaryDirectory() as tempdir:
    for file in input_files:
        shutil.copyfile(file, Path(tempdir) / Path(file).name)
    cmd_pars = [
        meta["executable"],
        "--input", tempdir,
        "--reference", reference,
        "--output", output,
        "---cpus", "2"
    ]
    subprocess.run(cmd_pars, check=True)

logger.info("> Check if file exists")

output_path = Path(output)
assert (output_path / "Log.final.out" ).is_file(), "No output log was created."
assert (output_path / "SJ.out.tab" ).is_file(), "No output file was created."


## Test 2: use input files
logger.info("> Running command with fastq files")
output = "test_output2"

cmd_pars = [
    meta["executable"],
    "--input", input_files[0],
    "--input", input_files[1],
    "--reference", reference,
    "--output", output,
    "---cpus", "8",
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")
assert path.exists(output + "/Log.final.out"), "No output log was created."
assert path.exists(output + "/SJ.out.tab"), "No output was created."


logger.info("> Completed Successfully!")