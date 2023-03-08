import subprocess
from os import path
import logging
from sys import stdout

## VIASH START
meta = {
    "functionality_name": "cellranger_mkfastq",
    "resources_dir": "resources_test"
}
## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("> Running command.")
input = meta["resources_dir"] + "/cellranger_tiny_bcl/bcl"
sample_sheet = meta["resources_dir"] + "/cellranger_tiny_bcl/bcl/sample_sheet.csv"
output = "test_output"

cmd_pars = [
    meta["executable"],
    "--input", input,
    "--sample_sheet", sample_sheet,
    "--output", output
]
if meta['cpus']:
    cmd_pars.extend(["---cpus", str(meta['cpus'])])
if meta['memory_gb']:
    cmd_pars.extend(["---memory", f"{meta['memory_gb']}GB"])
subprocess.check_call(cmd_pars, encoding="utf-8", timeout=500)

logger.info("> Check if file exists")
assert path.exists(output + "/H35KCBCXY/test_sample"), "No output was created."

logger.info("> Completed Successfully!")
