import subprocess
from os import path
import mudata as md
import numpy as np

input = meta["resources_dir"] + "/WTA.bd_rhapsody.output_raw"
output = "output1.h5mu"

cmd_pars = [
    meta["executable"],
    "--input", input,
    "--output", output,
    "--id", "foo",
    "--output_compression", "gzip",
]
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists(output), "No output was created."

# read it with scanpy
data = md.read_h5mu(output)

# check whether gex was found
assert np.array_equal(data.var["feature_types"].unique(), ["Gene Expression"]), "Output X should only contain Gene Expression vars."

# check whether gene was found
assert "PDE4DIP" in data.var_names, 'Output should contain gex column "PDE4DIP".'
