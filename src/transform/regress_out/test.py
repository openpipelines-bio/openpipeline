import sys
import pytest
import mudata as mu
import numpy as np

## VIASH START
meta = {"name": "lognorm", "resources_dir": "resources_test/"}
## VIASH END


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def input_data(input_path):
    return mu.read_h5mu(input_path)


@pytest.fixture
def input_h5mu(input_data):
    input_data.obs["var"] = np.random.rand(input_data.n_obs)
    input_data.mod["rna"].obs["var"] = input_data.obs["var"]
    input_data.mod["prot"].obs["var"] = input_data.obs["var"]
    input_data.mod["rna"].layers["input"] = input_data.mod["rna"].X
    return input_data


@pytest.fixture
def input_h5mu_path(write_mudata_to_file, input_h5mu):
    return write_mudata_to_file(input_h5mu)


@pytest.fixture
def output_h5mu_path(tmp_path):
    return tmp_path / "output.h5mu"


def test_regress_out(run_component, input_h5mu_path, output_h5mu_path):
    # execute command
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "var",
        "--output_compression",
        "gzip",
    ]
    run_component(cmd_pars)

    assert output_h5mu_path.is_file(), "No output was created."

    mu_input = mu.read_h5mu(input_h5mu_path)
    mu_output = mu.read_h5mu(output_h5mu_path)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["prot"].'
    assert "prot" in mu_output.mod, 'Output should contain data.mod["prot"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]
    prot_in = mu_input.mod["prot"]
    prot_out = mu_output.mod["prot"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert prot_in.shape == prot_out.shape, "Should have same shape as before"

    assert np.mean(rna_in.X) != np.mean(rna_out.X), "Expression should have changed"


def test_regress_out_with_layers(run_component, input_h5mu_path, output_h5mu_path):
    # execute command
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "var",
        "--input_layer",
        "input",
        "--output_layer",
        "output",
        "--output_compression",
        "gzip",
    ]
    run_component(cmd_pars)

    mu_input = mu.read_h5mu(input_h5mu_path)
    mu_output = mu.read_h5mu(output_h5mu_path)

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]

    assert np.mean(rna_in.layers["input"]) != np.mean(
        rna_out.layers["output"]
    ), "RNA expression should have changed"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
