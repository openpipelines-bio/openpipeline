from os import path
from mudata import read_h5mu
import pytest
import sys

## VIASH START
meta = {
    'executable': 'target/executable/correction/cellbender_remove_background/cellbender_remove_background',
    'resources_dir': 'resources_test/pbmc_1k_protein_v3'
}
## VIASH END

file_input = meta["resources_dir"] + "/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5mu"

def test_run(run_component, random_h5mu_path):
    print("> Check whether cellbender works when it should be working")
    
    file_output = random_h5mu_path()
    # run cellbender
    cmd_pars = [
        "--input", file_input,
        "--output", file_output,
        "--epochs", "5",
        "--output_compression", "gzip"
    ]
    # todo: if cuda is available, add --cuda
    run_component(cmd_pars)

    # check if file exists
    assert path.exists(file_output), "No output was created."

    data = read_h5mu(file_output)

    # check whether gex was found
    assert data.mod["rna"].var["feature_types"].unique() == [
        "Gene Expression"
    ], "Output X should only contain Gene Expression vars."

    # check whether ab counts were found
    assert "prot" in data.mod, 'Output should contain data.mod["rna"].'
  
  
if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))