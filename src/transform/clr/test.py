import sys
import pytest
from mudata import read_h5mu
import numpy as np

## VIASH START
meta = {
    'executable': 'target/executable/transform/clr/clr',
    'resources_dir': './resources_test/',
    'cpus': 2,
    'config': "./src/transform/clr/config.vsh.yaml"
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

def test_clr(run_component, tmp_path):
    output_file = tmp_path / "foo.h5mu"

    run_component([
        "--input", input_file,
        "--output", str(output_file),
        "--output_compression", "gzip",
        "--output_layer", "clr"
    ])
    assert output_file.is_file()
    output_h5mu = read_h5mu(output_file)
    assert 'clr' in output_h5mu.mod['prot'].layers.keys()
    assert output_h5mu.mod['prot'].layers['clr'] is not None
    input = read_h5mu(input_file)
    input_col = input.mod['prot'].X[:,0].toarray()
    result_col = output_h5mu.mod['prot'].layers['clr'][:,0].toarray()
    expected_col = np.log1p(input_col / np.exp(np.log1p(input_col).sum(axis=0) / input_col.size ))
    np.testing.assert_allclose(result_col, expected_col)


def test_clr_select_input_layer(run_component, tmp_path):
    output_file = tmp_path / "foo.h5mu"

    input_data = read_h5mu(input_file)
    input_data.mod['prot'].layers['test_layer'] = input_data.mod["prot"].X.copy()
    input_data.mod["prot"].X = None
    
    temp_input_file = tmp_path / "temp.h5mu"
    input_data.write(temp_input_file)

    run_component([
        "--input", temp_input_file,
        "--output", str(output_file),
        "--output_compression", "gzip",
        "--output_layer", "clr",
        "--input_layer", "test_layer",
    ])
    assert output_file.is_file()
    output_h5mu = read_h5mu(output_file)
    assert 'clr' in output_h5mu.mod['prot'].layers.keys()
    assert output_h5mu.mod['prot'].layers['clr'] is not None

def test_clr_output_to_x(run_component, tmp_path):
    output_file = tmp_path / "foo.h5mu"

    original_x = read_h5mu(input_file).mod['prot'].X
    run_component([
        "--input", input_file,
        "--output", str(output_file),
        "--output_compression", "gzip",
    ])
    assert output_file.is_file()
    output_h5mu = read_h5mu(output_file)
    assert 'clr' not in output_h5mu.mod['prot'].layers
    assert not np.all(np.isclose(original_x.toarray(), 
                                 output_h5mu.mod['prot'].X.toarray(), 
                                 rtol=1e-07, atol=1e-07))
    
def test_clr_set_axis(run_component, tmp_path):
    output_file = tmp_path / "foo.h5mu"

    run_component([
        "--input", input_file,
        "--output", str(output_file),
        "--output_compression", "gzip",
        "--output_layer", "clr",
        "--axis", "1",
    ])
    assert output_file.is_file()
    output_h5mu = read_h5mu(output_file)
    assert 'clr' in output_h5mu.mod['prot'].layers.keys()
    assert output_h5mu.mod['prot'].layers['clr'] is not None
    input = read_h5mu(input_file)
    input_row = input.mod['prot'].X[0].toarray()
    result_row = output_h5mu.mod['prot'].layers['clr'][0].toarray()
    expected_row = np.log1p(input_row / np.exp(np.log1p(input_row).sum(axis=1) / input_row.size ))
    np.testing.assert_allclose(result_row, expected_row)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))