import subprocess
import mudata as mu
import pytest
from pathlib import Path
from mudata import read_h5mu
import re
import numpy as np
import sys

## VIASH START
meta = {
    'functionality_name': 'foo',
    'resources_dir': 'resources_test/',
    'executable': './target/docker/dimred/pca/pca',
    'config': './src/dimred/pca/config.vsh.yaml'
}


## VIASH END

resources_dir, functionality_name = meta["resources_dir"], meta["functionality_name"]
input = meta["resources_dir"] + "pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"

def test_pca(run_component):
    run_component([
            "--input", input,
            "--output",  "output.h5mu",
            "--obsm_output", "X_foo",
            "--num_components", "26",
            "--overwrite"
        ])
    assert Path("output.h5mu").is_file()
    data = mu.read_h5mu("output.h5mu")

    # check whether pca was found
    assert data.mod["rna"].obsm["X_foo"].shape == (data.n_obs, 26)
    assert "highly_variable" not in data.mod["rna"].var.columns
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "pca_variance" in data.mod['rna'].uns
    assert "pca_loadings" in data.mod['rna'].varm
    assert "X_foo" in data.mod['rna'].obsm
    # GH298
    assert not np.array_equal(data.mod['rna'].uns['pca_variance']['variance'],
                              data.mod['rna'].uns['pca_variance']['variance_ratio'])

def test_no_overwrite_but_field_also_not_present(run_component, tmp_path):
    input_data = read_h5mu(input)
    input_data.mod['rna'].uns.pop('pca_variance')
    input_data.mod['rna'].varm.pop('pca_loadings')
    input_data.mod['rna'].obsm.pop('X_pca')
    tmp_file = tmp_path / "input_data_adjusted.h5mu"
    input_data.write(tmp_file)
    run_component([
        "--input", tmp_file,
        "--output",  "output.h5mu",
        "--obsm_output", "X_foo",
        "--num_components", "26",
        "--output_compression", "gzip"
    ])
    assert Path("output.h5mu").is_file()
    data = mu.read_h5mu("output.h5mu")

    # check whether pca was found
    assert data.mod["rna"].obsm["X_foo"].shape == (data.n_obs, 26)
    assert "highly_variable" not in data.mod["rna"].var.columns
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "pca_variance" in data.mod['rna'].uns
    assert "pca_loadings" in data.mod['rna'].varm
    assert "X_foo" in data.mod["rna"].obsm


def test_selecting_input_layer(run_component, tmp_path):
    input_data = read_h5mu(input)
    input_data.mod['rna'].layers['test'] = input_data.mod['rna'].X
    tmp_file = tmp_path / "input_data_adjusted.h5mu"
    input_data.write_h5mu(tmp_file)
    run_component([
            "--input",tmp_file,
            "--output",  "output.h5mu",
            "--obsm_output", "test_foo",
            "--num_components", "26",
            "--layer", "test",
            "--overwrite"
        ])
    assert Path("output.h5mu").is_file()
    data = mu.read_h5mu("output.h5mu")
    assert "test_foo" in data.mod["rna"].obsm
    assert data.mod["rna"].obsm["test_foo"].shape == (data.n_obs, 26)
    assert "highly_variable" not in data.mod["rna"].var.columns
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "pca_variance" in data.mod['rna'].uns
    assert "pca_loadings" in data.mod['rna'].varm
    assert "X_pca" in data.mod['rna'].obsm

def test_highly_variable_column_does_not_exist_raises(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
                "--input", input,
                "--output",  "output.h5mu",
                "--obsm_output", "X_foo",
                "--num_components", "26",
                "--var_input", "does_not_exist"
            ])
    error_found = re.search((r"ValueError: Requested to use \.var column does_not_exist as "
                             r"a selection of genes to run the PCA on, but the column is "
                             r"not available for modality rna"),
                             err.value.stdout.decode('utf-8'))
    assert error_found

def test_select_highly_variable_column(run_component, tmp_path):
    input_data = read_h5mu(input)
    input_data.mod['rna'].var["filter_with_hvg"] = True
    tmp_file = tmp_path / "input_data_adjusted.h5mu"
    input_data.write_h5mu(tmp_file)
    run_component([
            "--input", tmp_file,
            "--output",  "output.h5mu",
            "--obsm_output", "X_foo",
            "--num_components", "26",
            "--var_input", "filter_with_hvg",
            "--overwrite"
        ])
    assert Path("output.h5mu").is_file()
    data = mu.read_h5mu("output.h5mu")
    # check whether pca was found
    assert data.mod["rna"].obsm["X_foo"].shape == (data.n_obs, 26)
    assert "highly_variable" not in data.mod["rna"].var.columns
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "pca_variance" in data.mod['rna'].uns
    assert "pca_loadings" in data.mod['rna'].varm
    assert "X_pca" in data.mod['rna'].obsm

def test_raise_if_input_layer_is_missing(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
                "--input", input,
                "--output",  "output.h5mu",
                "--obsm_output", "X_foo",
                "--num_components", "26",
                "--layer", "does_not_exist",
                "--var_input", "filter_with_hvg"
            ])
    error_found = re.search(r'ValueError: does_not_exist was not found in modality rna\.',
                            err.value.stdout.decode('utf-8'))
    assert error_found

def test_output_field_already_present_raises(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input,
            "--output",  "output.h5mu",
            "--obsm_output", "X_foo",
            "--num_components", "26"
        ])
    error_found = re.search((r'ValueError: Requested to create field pca_loadings in \.varm for '
                             r'modality rna, but field already exists\.'),
                            err.value.stdout.decode('utf-8'))
    assert error_found

if __name__ == '__main__':
    sys.exit(pytest.main([__file__], plugins=["viashpy"]))