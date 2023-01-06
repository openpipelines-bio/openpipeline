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
@pytest.fixture
def viash_executable():
    return './bin/viash'
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

@pytest.fixture()
def lognormed_test_data(tmp_path, input_path):
    temp_h5mu = tmp_path / "lognormed.h5mu"
    mu_in = mu.read_h5mu(input_path)
    rna_in = mu_in.mod["rna"]
    assert "filter_with_hvg" not in rna_in.var.columns
    sc.pp.log1p(rna_in)
    mu_in.write_h5mu(temp_h5mu)
    return temp_h5mu


@pytest.fixture()
def filter_data(tmp_path, input_path):
    temp_h5mu = tmp_path / "filtered.h5mu"
    mu_in = mu.read_h5mu(input_path)
    rna_in = mu_in.mod["rna"]
    sc.pp.filter_genes(rna_in, min_counts=20)
    mu_in.write_h5mu(temp_h5mu)
    return temp_h5mu


def test_filter_with_hvg(run_component, lognormed_test_data):
    out = run_component([
        "--flavor", "seurat",
        "--input", lognormed_test_data,
        "--output", "output.h5mu"])
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

def test_filter_with_hvg_cell_ranger(run_component, filter_data):
    run_component([
        "--input", filter_data,
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
