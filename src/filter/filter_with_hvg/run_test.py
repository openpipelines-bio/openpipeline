import os
import subprocess
import scanpy as sc
import mudata as mu
import logging
import sys
import pytest
import re


## VIASH START
meta = {
    'resources_dir': 'resources_test/',
    'config': './src/filter/filter_with_hvg/config.vsh.yaml',
    'executable': './target/docker/filter/filter_with_hvg/filter_with_hvg'
}

## VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

@pytest.fixture
def input_data(input_path):
    mu_in = mu.read_h5mu(input_path)
    return mu_in

@pytest.fixture
def lognormed_test_data(input_data):
    rna_in = input_data.mod["rna"]
    assert "filter_with_hvg" not in rna_in.var.columns
    log_transformed = sc.pp.log1p(rna_in, copy=True)
    rna_in.layers['log_transformed'] = log_transformed.X
    rna_in.uns['log1p'] = log_transformed.uns['log1p']
    return input_data

@pytest.fixture
def lognormed_test_data_path(tmp_path, lognormed_test_data):
    temp_h5mu = tmp_path / "lognormed.h5mu"
    lognormed_test_data.write_h5mu(temp_h5mu)
    return temp_h5mu

@pytest.fixture
def lognormed_batch_test_data_path(tmp_path, lognormed_test_data):
    temp_h5mu = tmp_path / "lognormed_batch.h5mu"
    rna_mod = lognormed_test_data.mod['rna']
    rna_mod.obs['batch'] = 'A'
    column_index = rna_mod.obs.columns.get_indexer(['batch'])
    rna_mod.obs.iloc[slice(rna_mod.n_obs//2, None), column_index] = 'B'
    lognormed_test_data.write_h5mu(temp_h5mu)
    return temp_h5mu

@pytest.fixture()
def filter_data_path(tmp_path, input_data):
    temp_h5mu = tmp_path / "filtered.h5mu"
    rna_in = input_data.mod["rna"]
    sc.pp.filter_genes(rna_in, min_counts=20)
    input_data.write_h5mu(temp_h5mu)
    return temp_h5mu


def test_filter_with_hvg(run_component, lognormed_test_data_path):
    out = run_component([
        "--flavor", "seurat",
        "--input", lognormed_test_data_path,
        "--output", "output.h5mu",
        "--layer", "log_transformed"])
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns

def test_filter_with_hvg_batch_with_(run_component, lognormed_batch_test_data_path):
    """
    Make sure that selecting a layer works together with obs_batch_key.
    https://github.com/scverse/scanpy/issues/2396
    """
    run_component([
        "--flavor", "seurat",
        "--input", lognormed_batch_test_data_path,
        "--output", "output.h5mu",
        "--obs_batch_key", "batch",
        "--layer", "log_transformed"])
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns

def test_filter_with_hvg_seurat_v3_requires_n_top_genes(run_component, input_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_path,
            "--flavor", "seurat_v3", # Uses raw data.
            "--output", "output.h5mu"])
    assert re.search(f"When flavor is set to 'seurat_v3', you are required to set 'n_top_genes'\.",
                     err.value.stdout.decode('utf-8'))

def test_filter_with_hvg_seurat_v3(run_component, input_path):
    run_component([
        "--input", input_path,
        "--flavor", "seurat_v3", # Uses raw data.
        "--output", "output.h5mu",
        "--n_top_genes", "50"])
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns

def test_filter_with_hvg_cell_ranger(run_component, filter_data_path):
    run_component([
        "--input", filter_data_path,
        "--flavor", "cell_ranger", # Must use filtered data.
        "--output", "output.h5mu"])
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns

def test_filter_with_hvg_cell_ranger_unfiltered_data_change_error_message(run_component, input_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_path,
            "--flavor", "cell_ranger", # Must use filtered data, but in this test we use unfiltered data
            "--output", "output.h5mu"])
    assert re.search(r"Scanpy failed to calculate hvg\. The error "
                     r"returned by scanpy \(see above\) could be the "
                     r"result from trying to use this component on unfiltered data\.",
                    err.value.stdout.decode('utf-8'))


if __name__ == '__main__':
    sys.exit(pytest.main([__file__], plugins=["viashpy"]))
