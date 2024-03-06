import sys
import pytest
from uuid import uuid4
import numpy as np
from mudata import read_h5mu

## VIASH START
meta = {
    "executable": "./target/docker/scaling/scale",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
    "config": "./src/scaling/config.vsh.yaml",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

@pytest.fixture
def random_path(tmp_path):
    def wrapper(extension=None):
        extension = "" if not extension else extension
        return tmp_path / f"{uuid4()}.{extension}"
    return wrapper 


def test_scaling(run_component, random_path):
    """
    Output data must be centered around mean 0 and it has unit variance.
    """
    output_file = random_path()

    run_component([
        "--input", input_file,
        "--output", output_file,
        "--ouput_compression", "gzip"])

    assert output_file.is_file()
    output_data = read_h5mu(output_file)
    output_x = output_data.mod['rna'].X
    mean = np.mean(output_x, axis=0, dtype=np.float64)
    variance = np.multiply(output_x, output_x).mean(axis=0, dtype=np.float64) - mean**2
    variance[variance == 0] = 1
    assert np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07))
    assert np.all(np.isclose(variance, 1, rtol=1e-03, atol=1e-03))

def test_scaling_input_layer(run_component, random_path):
    output_file = random_path()

    run_component([
        "--input", input_file,
        "--output", output_file,
        "--input_layer", 
        "--ouput_compression", "gzip"]) 

def test_scaling_noncenter(run_component, random_path):
    """
    Check if centering can be disabled.
    """
    output_file = random_path()

    run_component([
        "--input", input_file,
        "--output", str(output_file),
        "--zero_center", "false"])
    assert output_file.is_file()
    output_data = read_h5mu(output_file)
    output_x = output_data.mod['rna'].X
    mean = np.mean(output_x, axis=0, dtype=np.float64)
    assert not np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07))

def test_scaling_maxvalue(run_component, random_path):
    """
    Check if output data is clipped when using --max_value
    """
    output_file = random_path()

    run_component([
        "--input", input_file,
        "--output", str(output_file),
        "--max_value", "0.5"])
    assert output_file.is_file()
    output_data = read_h5mu(output_file)
    output_x = output_data.mod['rna'].X
    assert np.all(output_x <= 0.5)

def test_scaling_modality(run_component, random_path):
    """
    Check if 'rna' modality remain untouched when using '--modality prot' argument.
    """
    output_file = random_path()

    run_component([
        "--input", input_file,
        "--output", str(output_file),
        "--modality", "prot"])
    assert output_file.is_file()
    input_data =  read_h5mu(input_file)
    output_data = read_h5mu(output_file)
    output_rna = output_data.mod['rna'].X
    assert np.allclose(input_data.mod['rna'].X.todense(), output_rna.todense(), equal_nan=True)

    output_prot =  output_data.mod['prot'].X
    mean = np.mean(output_prot, axis=0, dtype=np.float64)
    variance = np.multiply(output_prot, output_prot).mean(axis=0, dtype=np.float64) - mean**2
    variance[variance == 0] = 1
    assert np.all(np.isclose(mean, 0, rtol=1e-07, atol=1e-07))
    assert np.all(np.isclose(variance, 1, rtol=1e-03, atol=1e-03))

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
