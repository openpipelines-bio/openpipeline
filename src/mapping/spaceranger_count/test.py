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
    "functionality_name": "spaceranger_count",
    "resources_dir": "resources_test"
}
## VIASH END

logger.info("> Running command with folder")
input = meta["resources_dir"] + "/spaceranger_tiny_fastq/fastq/"
reference = meta["resources_dir"] + "/spaceranger_tiny_fastq/reference/refdata-gex-mm10-2020-A"
image = meta["resources_dir"] + "/spaceranger_tiny_fastq/image/V1_Adult_Mouse_Brain_image.tif"
output = "test_output"

cmd_pars = [
    meta["executable"],
    "--input", input,
    "--reference", reference,
    "--image", image,
    "--slide", "V19L01-041",
    "--area", "C1",
    "--output", output,
    "---cores", "2",
    "---memory", "5gb"
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

logger.info("> Check if file exists")
assert path.exists(output + "/filtered_feature_bc_matrix.h5"), "No output was created."

logger.info("> Completed Successfully!")