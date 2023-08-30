import pytest
import sys
from mudata import read_h5mu
from pathlib import Path
from subprocess import CalledProcessError
import re

## VIASH START
meta = {
    'functionality_name': './target/native/transform/delete_layer/delete_layer',
    'resources_dir': './resources_test/'
}
## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input_file = f"{resources_dir}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


def test_delete_layer(run_component, tmp_path):
    original_input_data = read_h5mu(input_file)
    new_layer = original_input_data.mod['rna'].X
    original_input_data.mod['rna'].layers['test'] = new_layer
    tempfile = tmp_path / "input.h5mu"
    original_input_data.write_h5mu(tempfile.name)
    run_component([
        "--input", tempfile.name,
        "--modality", "rna",
        "--layer", "test",
        "--output", "deleted_layer.h5mu"])
    assert Path("deleted_layer.h5mu").is_file()
    output_data = read_h5mu('deleted_layer.h5mu')
    assert 'test' not in output_data.mod['rna'].layers.keys()
    assert set(output_data.mod) == {'rna', 'prot'}

def test_missing_layer_raises(run_component):
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--modality", "rna",
            "--layer", "test",
            "--output", "missing_layer.h5mu"])
    re.search(r"Layer 'test' is not present in modality rna\.",
              err.value.stdout.decode('utf-8'))

def test_missing_layer_missing_ok(run_component):
    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--layer", "test",
        "--output", "missing_layer_ok.h5mu",
        "--missing_ok"])
    assert Path("missing_layer_ok.h5mu").is_file()
    output_data = read_h5mu('missing_layer_ok.h5mu')
    assert 'test' not in output_data.mod['rna'].layers.keys()
    assert set(output_data.mod) == {'rna', 'prot'}

@pytest.mark.parametrize("output_compression", ["gzip", "lzf"])
def test_delete_layer_with_compression(run_component, tmp_path, output_compression):
    original_input_data = read_h5mu(input_file)
    new_layer = original_input_data.mod['rna'].X
    original_input_data.mod['rna'].layers['test'] = new_layer
    tempfile = tmp_path / "input.h5mu"
    original_input_data.write_h5mu(tempfile.name)
    run_component([
        "--input", tempfile.name,
        "--modality", "rna",
        "--layer", "test",
        "--output", "deleted_layer.h5mu",
        "--output_compression", output_compression])
    assert Path("deleted_layer.h5mu").is_file()
    output_data = read_h5mu('deleted_layer.h5mu')
    assert 'test' not in output_data.mod['rna'].layers.keys()
    assert set(output_data.mod) == {'rna', 'prot'}


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))