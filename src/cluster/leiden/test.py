import mudata as mu
import pytest
import sys

## VIASH START
meta = {
    "name": "foo",
    "resources_dir": "resources_test/",
    "cpus": 2,
    "config": "./src/cluster/leiden/config.vsh.yaml",
    "executable": "./target/executable/cluster/leiden/leiden",
}
## VIASH END


@pytest.fixture()
def input_path():
    return meta["resources_dir"] + "/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture()
def input_data(input_path):
    return mu.read_h5mu(input_path)


@pytest.fixture()
def mudata_custom_connectivities_key(input_data, random_h5mu_path):
    result = input_data.copy()
    result.mod["rna"].obsp["custom_connectivities"] = (
        result.mod["rna"].obsp["connectivities"].copy()
    )
    del result.mod["rna"].obsp["connectivities"]
    output_path = random_h5mu_path()
    result.write(output_path)
    return output_path


# @pytest.fixture
# def random_h5mu_path(tmp_path):
#     def wrapper():
#         unique_filename = f"{str(uuid.uuid4())}.h5mu"
#         temp_file = tmp_path / unique_filename
#         return temp_file
#     return wrapper


@pytest.mark.parametrize("compression", ["gzip", ""])
@pytest.mark.parametrize(
    "output_key,expected_output_key", [("fooleiden", "fooleiden"), ("", "leiden")]
)
def test_leiden(
    input_path,
    run_component,
    random_h5mu_path,
    compression,
    output_key,
    expected_output_key,
):
    output_path = random_h5mu_path()
    args = ["--input", input_path, "--resolution", "1;0.25", "--output", output_path]
    if compression:
        args.extend(["--output_compression", compression])
    if output_key:
        args.extend(["--obsm_name", output_key])

    run_component(args)
    assert output_path.exists(), "No output was created."
    data = mu.read_h5mu(output_path)
    assert expected_output_key in data.mod["rna"].obsm, (
        f"Expected to find key '{expected_output_key}' in .obsm"
    )
    # check whether leiden.custom.resolution was found
    assert "1.0" in data.mod["rna"].obsm[expected_output_key].columns, (
        "Output should contain resolution 1.0."
    )
    assert "0.25" in data.mod["rna"].obsm[expected_output_key].columns, (
        "Output should contain resolution 0.25."
    )


def test_leiden_custom_connectivities_key(
    mudata_custom_connectivities_key, run_component, random_h5mu_path
):
    output_path = random_h5mu_path()
    run_component(
        [
            "--input",
            mudata_custom_connectivities_key,
            "--obsm_name",
            "fooleiden",
            "--resolution",
            "1;0.25",
            "--output",
            output_path,
            "--obsp_connectivities",
            "custom_connectivities",
            "--output_compression",
            "gzip",
        ]
    )
    assert output_path.exists(), "No output was created."
    data = mu.read_h5mu(output_path)
    # check whether leiden.custom.resolution was found
    assert "1.0" in data.mod["rna"].obsm["fooleiden"].columns, (
        "Output should contain resolution 1.0."
    )
    assert "0.25" in data.mod["rna"].obsm["fooleiden"].columns, (
        "Output should contain resolution 0.25."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
