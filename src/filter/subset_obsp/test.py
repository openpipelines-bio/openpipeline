import sys
import pytest
import mudata as mu

## VIASH START
meta = {"resources_dir": "resources_test/pbmc_1k_protein_v3/"}
## VIASH END


@pytest.fixture
def input_h5mu():
    input = mu.read_h5mu(f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu")
    input.mod["rna"].obs["filter_column"] = "group_2"
    input.mod["rna"].obs["filter_column"][:50] = "group_1"
    return input


@pytest.fixture
def input_path(write_mudata_to_file, input_h5mu):
    return write_mudata_to_file(input_h5mu)


def test_subset_obsp(input_path, run_component, tmp_path):
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component(
        [
            "--input",
            input_path,
            "--output",
            str(output_path),
            "--input_obsp_key",
            "distances",
            "--input_obs_key",
            "filter_column",
            "--input_obs_value",
            "group_1",
            "--output_obsm_key",
            "group_1",
        ]
    )

    assert output_path.is_file(), "Output file not found"

    # check output file
    mu_out = mu.read_h5mu(output_path)

    assert "group_1" in mu_out.mod["rna"].obsm, "Output should contain group_1 in .obsm"
    assert mu_out.mod["rna"].obsm["group_1"].shape[1] == 50, (
        "Obsm should only contain a subset of the original obsp matrix"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
