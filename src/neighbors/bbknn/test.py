import sys
import pytest
from mudata import read_h5mu

## VIASH START
meta = {
    "executable": "./target/executable/neighbors/bbknn/bbknn",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def sample_mudata(tmp_path):
    tmp_input_path = tmp_path / "input.h5mu"

    # create input data
    mudata = read_h5mu(input_file)
    rna_adata = mudata.mod["rna"]

    # remove previous output (if any)
    if "connectivities" in rna_adata.obsp:
        del rna_adata.obsp["connectivities"]
    if "distances" in rna_adata.obsp:
        del rna_adata.obsp["distances"]
    if "neighbors" in rna_adata.uns:
        del rna_adata.uns["neighbors"]

    # write to file
    mudata.write(tmp_input_path)

    return tmp_input_path, mudata


def test_simple_integration(run_component, tmp_path, sample_mudata):
    tmp_input_path, mudata = sample_mudata
    output_path = tmp_path / "output.h5mu"
    print(mudata, flush=True)
    # run component
    run_component(
        [
            "--input",
            str(tmp_input_path),
            "--output",
            str(output_path),
            "--obs_batch",
            "harmony_integration_leiden_1.0",
            "--obsm_input",
            "X_pca",
            "--output_compression",
            "gzip",
        ]
    )
    assert output_path.exists()
    data = read_h5mu(output_path).mod["rna"]
    assert "connectivities" in data.obsp
    assert "distances" in data.obsp
    assert "neighbors" in data.uns


def test_alternative_names(run_component, tmp_path, sample_mudata):
    tmp_input_path, mudata = sample_mudata
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component(
        [
            "--input",
            str(tmp_input_path),
            "--output",
            str(output_path),
            "--obs_batch",
            "harmony_integration_leiden_1.0",
            "--obsm_input",
            "X_pca",
            "--output_compression",
            "gzip",
            "--uns_output",
            "my_neighbors",
            "--obsp_connectivities",
            "my_connectivities",
            "--obsp_distances",
            "my_distances",
        ]
    )
    assert output_path.exists()
    data = read_h5mu(output_path).mod["rna"]
    assert "my_connectivities" in data.obsp
    assert "my_distances" in data.obsp
    assert "my_neighbors" in data.uns


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
