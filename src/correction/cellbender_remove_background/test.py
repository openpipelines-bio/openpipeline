from os import path, listdir
from mudata import read_h5mu
import subprocess

## VIASH START
meta = {
    "executable": "target/executable/correction/cellbender_remove_background/cellbender_remove_background",
    "resources_dir": "resources_test/pbmc_1k_protein_v3",
}
## VIASH END

file_input = (
    meta["resources_dir"] + "/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)
file_output = "output.h5mu"
dir_output_raw = "cellbender_output"

print("> Check whether cellbender works when it should be working")

# run cellbender
cmd_pars = [
    meta["executable"],
    "--input",
    file_input,
    "--output",
    file_output,
    "--output_raw",
    dir_output_raw,
    "--epochs",
    "5",
    "--output_compression",
    "gzip",
]
# todo: if cuda is available, add --cuda
out = subprocess.check_output(cmd_pars).decode("utf-8")

# check if file exists
assert path.exists(file_output), "No output was created."

# check whether the full cellbender output bundle was published
assert path.isdir(dir_output_raw), "No output bundle was created."
assert listdir(dir_output_raw), "Output bundle is empty."
expected_bundle_files = [
    "output_cell_barcodes.csv",
    "output_filtered.h5",
    "output_metrics.csv",
    "output_posterior.h5",
    "output_report.html",
    "output.h5",
    "output.pdf",
]
for bundle_file in expected_bundle_files:
    assert path.exists(path.join(dir_output_raw, bundle_file)), (
        f"Output bundle should contain '{bundle_file}'."
    )

data = read_h5mu(file_output)

# check whether gex was found
assert data.mod["rna"].var["feature_types"].unique() == ["Gene Expression"], (
    "Output X should only contain Gene Expression vars."
)

# check whether ab counts were found
assert "prot" in data.mod, 'Output should contain data.mod["rna"].'
