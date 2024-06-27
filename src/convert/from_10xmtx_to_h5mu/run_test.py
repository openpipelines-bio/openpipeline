from os import path
from mudata import read_h5mu
import pytest
import sys

## VIASH START
meta = {
    'resources_dir': 'resources_test/',
    'config': './src/convert/from_10xmtx_to_h5mu/config.vsh.yaml',
    'executable': './target/executable/convert/from_10xmtx_to_h5mu/from_10xmtx_to_h5mu',
}
## VIASH END

input = meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix"

def test_run(run_component, random_h5mu_path):
    output = random_h5mu_path()
    cmd_pars = [
        "--input", input,
        "--output", output,
        "--output_compression", "gzip",
    ]
    run_component(cmd_pars)

    # check if file exists
    assert path.exists(output), "No output was created."

    # read it with scanpy
    data = read_h5mu(output)

    # check whether gex was found
    assert data.mod["rna"].var["feature_types"].unique() == [
        "Gene Expression"
    ], "Output X should only contain Gene Expression vars."

    # check whether ab counts were found
    assert "prot" in data.mod, 'Output should contain data.mod["rna"].'

    # check whether gene was found
    assert (
        "CD3" in data.mod["prot"].var_names
    ), 'Output should contain antibody column "CD3".'
    
    
if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
