import os
import subprocess
import scanpy as sc
import muon

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/'
}
## VIASH END

orig_input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
input_path = "lognormed.h5mu"
output_path = "output.h5mu"

print("> Prepare test data")
mu_in = muon.read_h5mu(orig_input_path)
rna_in = mu_in.mod["rna"]
assert "filter_with_hvg" not in rna_in.var.columns
sc.pp.log1p(rna_in)
mu_in.write_h5mu(input_path)


print("> Run component")
out = subprocess.check_output(
    [
        f"./{meta['functionality_name']}", 
        "--input", input_path,
        "--output", output_path
    ]
).decode("utf-8")

print("> Check output file exists")
assert os.path.exists(output_path)
data = muon.read_h5mu(output_path)

print("> Check whether column has been added")
assert "filter_with_hvg" in data.mod["rna"].var.columns