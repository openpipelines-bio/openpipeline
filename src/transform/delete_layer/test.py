from subprocess import CalledProcessError
import sys
import pytest
from mudata import read_h5mu
from pathlib import Path

## VIASH START
meta = {
    'functionality_name': './target/native/transform/delete_layer/delete_layer',
    'resources_dir': './resources_test/'
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

def test_delete_layer(run_component, tmp_path):
    original_input_data = read_h5mu(input_file)
    new_layer = original_input_data.mod['rna'].X
    original_input_data.mod['rna'].layers['test'] = new_layer

    tempfile = tmp_path / "temp.h5mu"
    original_input_data.write_h5mu(tempfile.name)
    run_component([
        "--input", tempfile.name,
        "--modality", "rna",
        "--layer", "test",
        "--output", "deleted_layer.h5mu",
        "--output_compression", "gzip"])
    
    assert Path("deleted_layer.h5mu").is_file()
    output_data = read_h5mu('deleted_layer.h5mu')
    assert 'test' not in output_data.mod['rna'].layers.keys()

def test_missing_layer_raises(run_component, tmp_path):
    tmpfile = tmp_path / "temp.h5mu"
    with pytest.raises(CalledProcessError, match=r"Layer 'test' is not present in modality rna\."):
        run_component([
            "--input", input_file,
            "--modality", "rna",
            "--layer", str(tmpfile),
            "--output", "missing_layer.h5mu"],
            expected_raise=True)
    assert not tmpfile.is_file()

def test_missing_layer_missing_ok(run_component, tmp_path):
    tmpfile = tmp_path / "temp.h5mu"
    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--layer", "test",
        "--output", str(tmpfile),
        "--missing_ok"])
    assert tmpfile.is_file()
    output_data = read_h5mu(tmpfile)
    assert 'test' not in output_data.mod['rna'].layers.keys()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))