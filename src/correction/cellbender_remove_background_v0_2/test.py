from os import path
import muon as mu
import pytest

## VIASH START
meta = {
    "executable": "target/executable/correction/cellbender_remove_background/cellbender_remove_background",
    "resources_dir": "resources_test/pbmc_1k_protein_v3",
}
## VIASH END

file_raw = meta["resources_dir"] + "/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5"


@pytest.fixture
def subsampled_input(write_mudata_to_file):
    mdat = mu.read_10x_h5(file_raw)
    mdat = mdat[0:100000,]
    return write_mudata_to_file(mdat)


def test_run(run_component, random_h5mu_path, subsampled_input):
    print("> Check whether cellbender works when it should be working")

    # run cellbender
    output_file = random_h5mu_path()
    cmd_pars = [
        "--input",
        subsampled_input,
        "--output",
        output_file,
        "--epochs",
        "5",
        "--output_compression",
        "gzip",
    ]
    # todo: if cuda is available, add --cuda
    run_component(cmd_pars)

    # check if file exists
    assert path.exists(output_file), "No output was created."

    # read it with scanpy
    data = mu.read_h5mu(output_file)

    # check whether gex was found
    assert data.mod["rna"].var["feature_types"].unique() == [
        "Gene Expression"
    ], "Output X should only contain Gene Expression vars."

    # check whether ab counts were found
    assert "prot" in data.mod, 'Output should contain data.mod["rna"].'
