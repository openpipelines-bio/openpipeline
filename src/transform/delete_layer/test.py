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
    tempfile = tmp_path / "input.h5mu"
    output = tmp_path / "deleted_layer.h5mu"

    # create input file
    input_data = read_h5mu(input_file)
    new_layer = input_data.mod['rna'].X
    input_data.mod['rna'].layers['test'] = new_layer
    assert "test" in input_data.mod['rna'].layers.keys()
    input_data.write_h5mu(tempfile.name)

    # run command
    run_component([
        "--input", str(tempfile),
        "--modality", "rna",
        "--layer", "test",
        "--output", str(output)])
    
    # check if output is correct
    assert output.is_file()
    output_data = read_h5mu(output)
    assert 'test' not in output_data.mod['rna'].layers.keys()
    assert set(output_data.mod) == {'rna', 'prot'}

def test_missing_layer_raises(run_component, tmp_path):
    tmpfile = tmp_path / "temp.h5mu"
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--modality", "rna",
            "--layer", "test",
            "--output", str(tmpfile)])
    assert not tmpfile.is_file()
    assert re.search(r"Layer 'test' is not present in modality rna\.",
              err.value.stdout.decode('utf-8'))

def test_missing_layer_missing_ok(run_component, tmp_path):
    tmp_file = tmp_path / "temp.h5mu"
    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--layer", "test",
        "--output", str(tmp_file),
        "--missing_ok"])
    assert tmp_file.is_file()
    output_data = read_h5mu(tmp_file)
    assert 'test' not in output_data.mod['rna'].layers.keys()
    assert set(output_data.mod) == {'rna', 'prot'}

@pytest.mark.parametrize("output_compression", ["gzip", "lzf"])
def test_delete_layer_with_compression(run_component, tmp_path, output_compression):
    original_input_data = read_h5mu(input_file)
    new_layer = original_input_data.mod['rna'].X
    original_input_data.mod['rna'].layers['test'] = new_layer
    tempfile = tmp_path / "input.h5mu"
    original_input_data.write_h5mu(tempfile.name)
    output = tmp_path / "deleted_layer.h5mu"
    run_component([
        "--input", tempfile.name,
        "--modality", "rna",
        "--layer", "test",
        "--output", str(output),
        "--output_compression", output_compression])
    assert output.is_file()
    output_data = read_h5mu(output)
    assert 'test' not in output_data.mod['rna'].layers.keys()
    assert set(output_data.mod) == {'rna', 'prot'}


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))