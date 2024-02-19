import sys
import pytest
import re

from mudata import read_h5mu
from subprocess import CalledProcessError

## VIASH START
meta = {
    'functionality_name': './target/native/transform/delete_layer/delete_layer',
    'resources_dir': './resources_test/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

def test_delete_layer(run_component, tmp_path):
    temp_input = tmp_path / "input.h5mu"
    temp_output = tmp_path / "output.h5mu"

    # create input file
    input = read_h5mu(input_file)
    new_layer = input.mod['rna'].X
    input.mod['rna'].layers['test'] = new_layer
    assert "test" in input.mod['rna'].layers.keys()
    input.write_h5mu(temp_input)

    # run command
    run_component([
        "--input", str(temp_input),
        "--modality", "rna",
        "--layer", "test",
        "--output", str(temp_output)])
    
    # check if output is correct
    assert temp_output.is_file()
    output = read_h5mu(temp_output)
    assert 'test' not in output.mod['rna'].layers.keys()
    assert set(output.mod) == {'rna', 'prot'}

def test_missing_layer_raises(run_component, tmp_path):
    output = tmp_path / "temp.h5mu"
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--modality", "rna",
            "--layer", "test",
            "--output", str(output)])
    assert not output.is_file()
    assert "Layer 'test' is not present in modality rna." in \
            err.value.stdout.decode('utf-8')

def test_missing_layer_missing_ok(run_component, tmp_path):
    output = tmp_path / "temp.h5mu"
    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--layer", "test",
        "--output", str(output),
        "--missing_ok"])
    assert output.is_file()
    output_data = read_h5mu(output)
    assert 'test' not in output_data.mod['rna'].layers.keys()
    assert set(output_data.mod) == {'rna', 'prot'}

@pytest.mark.parametrize("output_compression", ["gzip", "lzf"])
def test_delete_layer_with_compression(run_component, tmp_path, output_compression):
    temp_input = tmp_path / "input.h5mu"
    output = tmp_path / "output.h5mu"

    # create temp input with 'test' layer
    original_input_data = read_h5mu(input_file)
    new_layer = original_input_data.mod['rna'].X
    original_input_data.mod['rna'].layers['test'] = new_layer
    original_input_data.write_h5mu(temp_input)

    # run component
    run_component([
        "--input", str(temp_input),
        "--modality", "rna",
        "--layer", "test",
        "--output", str(output),
        "--output_compression", output_compression])
    
    # check if output is correct
    assert output.is_file()
    output_data = read_h5mu(output)
    assert 'test' not in output_data.mod['rna'].layers.keys()
    assert set(output_data.mod) == {'rna', 'prot'}


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))