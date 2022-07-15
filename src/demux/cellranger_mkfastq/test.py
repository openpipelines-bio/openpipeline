import subprocess
from os import path

## VIASH START
meta = {
    "functionality_name": "cellranger_mkfastq",
    "resources_dir": "resources_test"
}
## VIASH END

print("> Running command", flush=True)
input = meta["resources_dir"] + "/cellranger_tiny_bcl/bcl"
sample_sheet = meta["resources_dir"] + "/cellranger_tiny_bcl/bcl/sample_sheet.csv"
output = "test_output"

cmd_pars = [
    "./" + meta["functionality_name"],
    "--input", input,
    "--sample_sheet", sample_sheet,
    "--output", output,
    "--memory", "5" 
]
subprocess.check_call(cmd_pars, encoding="utf-8", timeout=300)
print("> Check if file exists")
assert path.exists(output + "/H35KCBCXY/test_sample"), "No output was created."

print("> Completed Successfully!")