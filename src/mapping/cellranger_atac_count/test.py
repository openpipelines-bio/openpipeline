import subprocess
from os import path
import sys
from itertools import zip_longest, chain

## VIASH START
meta = {
    "name": "cellranger_atac_count",
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

logger.info("> Running command with folder")
input = meta["resources_dir"] + "/cellranger_atac_tiny_bcl/fastqs/HJN3KBCX2/test_sample/"
reference = meta["resources_dir"] + "/reference_gencodev41_chr1/reference_cellranger.tar.gz"
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
try:
    out = subprocess.check_output(cmd_pars).decode("utf-8")
except subprocess.CalledProcessError as e:
    logger.error(e.output)
    raise e
logger.info("> Check if file exists")
assert path.exists(output + "/filtered_peak_bc_matrix.h5"), "No output was created."
assert path.exists(output + "/fragments.tsv.gz"), "No fragments file was created."


logger.info("> Running command with fastq files")
# test_sample_S1_L001_R2_001.fastq.gz test_sample_S1_L001_R1_001.fastq.gz  test_sample_S1_L001_R3_001.fastq.gz
input_files = [
    input + "test_sample_S1_L001_I1_001.fastq.gz",
    input + "test_sample_S1_L001_R1_001.fastq.gz",
    input + "test_sample_S1_L001_R2_001.fastq.gz",
    input + "test_sample_S1_L001_R3_001.fastq.gz",    
]
output = "test_output2"


cmd_pars = [
    meta["executable"],
    *chain.from_iterable([("--input", input_file) for input_file in input_files]),
    "--reference", reference,
    "--output", output
]
if meta.get("cpus"):
    cmd_pars.extend(["---cpus", str(meta["cpus"])])
if meta.get("memory_gb"):
    cmd_pars.extend(["---memory", f"{meta['memory_gb']}gb"])
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")
assert path.exists(output + "/filtered_peak_bc_matrix.h5"), "No output was created."
assert path.exists(output + "/fragments.tsv.gz"), "No fragments file was created."

logger.info("> Completed Successfully!")