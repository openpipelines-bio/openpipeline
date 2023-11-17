import sys
import pytest
from mudata import read_h5mu
import numpy as np

## VIASH START
meta = {
    'executable': './target/docker/integrate/scpoli/scpoli',
    'resources_dir': './resources_test/pbmc_1k_protein_v3/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

@pytest.fixture
def input_with_batch(tmp_path):
    tmp_input_path = tmp_path / "input.h5mu"

    input_data = read_h5mu(input_file)

    mod = input_data.mod['rna']

    mod.obs["batch_1"] = np.random.choice(["A1", "A2"], mod.n_obs)
    mod.obs["batch_2"] = np.random.choice(["B1", "B2"], mod.n_obs)

    mod.obs["cell_type"] = np.random.choice(["celltype_A", "celltype_B"], mod.n_obs)

    input_data.write(tmp_input_path)

    return tmp_input_path, input_data

def test_simple_integration(run_component, input_with_batch, tmp_path):
    tmp_input_path, _ = input_with_batch
    output_path = tmp_path / "output.h5mu"

    # run component
    run_component([
        "--input", str(tmp_input_path),
        "--output", str(output_path),
        "--condition_keys", "batch_1",
        "--condition_keys", "batch_2",
        "--cell_type_keys", "cell_type",
        "--output_compression", "gzip", 
        "--n_epochs", "1"])
    assert output_path.is_file()

    # check output
    data = read_h5mu(output_path)
    assert "X_scpoli_integrated" in data.mod['rna'].obsm

def test_output_changes(run_component, input_with_batch, tmp_path):
    tmp_input_path, _ = input_with_batch
    output_path = tmp_path / "output.h5mu"
    model_dir_output_path = tmp_path / "output_model/"

    # run component
    run_component([
        "--input", str(tmp_input_path),
        "--output", str(output_path),
        "--condition_keys", "sample_id",
        "--cell_type_keys", "sample_id", 
        "--n_epochs", "1",
        "--obsm_output", "X_test", 
        "--output_model", str(model_dir_output_path)
        ])
    assert output_path.is_file()

    # check output
    data = read_h5mu(output_path)
    #assert "X_test" in data.mod['rna'].obsm
    assert "X_test" in data.mod['rna'].obsm


    assert model_dir_output_path.is_dir()

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))