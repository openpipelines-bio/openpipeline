import mudata as mu
import sys
from pathlib import Path
import pytest
import numpy as np
from subprocess import CalledProcessError

## VIASH START
meta = {
    'executable': './target/docker/filter/delimit_fraction/delimit_fraction',
    'resources_dir': 'resources_test/',
    'config': "./src/filter/delimit_fraction/config.vsh.yaml"
}

## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()


@pytest.fixture
def original_input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

@pytest.fixture
def input_h5mu(original_input_path):
    input_data = mu.read_h5mu(original_input_path)
    input_data.mod['rna'].obs['test_fraction'] = np.random.rand(input_data.mod['rna'].n_obs)
    return input_data


@pytest.fixture
def input_h5mu_string_data(original_input_path):
    input_data = mu.read_h5mu(original_input_path)
    input_data.mod['rna'].obs['test_fraction'] = np.random.choice(['these', 'are', 'random', 'values'], input_data.mod['rna'].n_obs)
    return input_data

@pytest.fixture
def input_path(input_h5mu, tmp_path):
    output_path = tmp_path / "temp_h5mu.h5mu"
    input_h5mu.write(output_path)
    return output_path


@pytest.fixture
def input_path_string_data(input_h5mu_string_data, tmp_path):
    output_path = tmp_path / "temp_h5mu_string.h5mu"
    input_h5mu_string_data.write(output_path)
    return output_path

@pytest.fixture
def input_n_rna_obs(input_h5mu):
    return input_h5mu.mod['rna'].n_obs

@pytest.fixture
def input_n_prot_obs(input_h5mu):
    return input_h5mu.mod['prot'].n_obs

@pytest.fixture
def input_n_rna_vars(input_h5mu):
    return input_h5mu.mod['rna'].n_vars

@pytest.fixture
def input_n_prot_vars(input_h5mu):
    return input_h5mu.mod['prot'].n_vars

def test_filter_nothing(run_component, input_path,
                        input_n_rna_obs, input_n_prot_obs,
                        input_n_rna_vars, input_n_prot_vars):
    run_component([
        "--input", input_path,
        "--output", "output-1.h5mu",
        "--min_fraction", "0",
        "--max_fraction", "1",
        "--output_compression", "gzip",
        "--obs_name_filter", "test_output",
        "--obs_fraction_column", "test_fraction"
        ])
    assert Path("output-1.h5mu").is_file()
    mu_out = mu.read_h5mu("output-1.h5mu")
    assert "test_output" in mu_out.mod["rna"].obs
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs == input_n_rna_obs
    assert new_vars == input_n_rna_vars
    assert mu_out.mod['prot'].n_obs == input_n_prot_obs
    assert mu_out.mod['prot'].n_vars == input_n_prot_vars
    assert mu_out.mod['rna'].obs['test_output'].all()
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"]
    assert list(mu_out.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"]

def test_filtering_a_little(run_component, input_path,
                            input_n_rna_obs, input_n_prot_obs,
                            input_n_rna_vars, input_n_prot_vars):
    run_component([
        "--input", input_path,
        "--output", "output-2.h5mu",
        "--min_fraction", "0.5",
        "--max_fraction", "0.7",
        "--output_compression", "gzip",
        "--obs_name_filter", "test_output",
        "--obs_fraction_column", "test_fraction"     
    ])

    assert Path("output-2.h5mu").is_file()
    mu_out = mu.read_h5mu("output-2.h5mu")
    new_obs = mu_out.mod['rna'].n_obs
    new_vars = mu_out.mod['rna'].n_vars
    assert new_obs == input_n_rna_obs
    assert new_vars == input_n_rna_vars
    assert mu_out.mod['prot'].n_obs == input_n_prot_obs
    assert mu_out.mod['prot'].n_vars == input_n_prot_vars
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"]
    assert not mu_out.mod['rna'].obs['test_output'].all()
    assert list(mu_out.mod['rna'].var['feature_types'].cat.categories) == ["Gene Expression"]
    assert list(mu_out.mod['prot'].var['feature_types'].cat.categories) == ["Antibody Capture"]


def test_filtering_wrong_data_raises(run_component, input_path_string_data):
    with pytest.raises(CalledProcessError) as err:
        run_component([
            "--input", input_path_string_data,
            "--output", "output-2.h5mu",
            "--min_fraction", "0.5",
            "--max_fraction", "0.7",
            "--output_compression", "gzip",
            "--obs_name_filter", "test_output",
            "--obs_fraction_column", "test_fraction"     
        ])
    assert "Column 'test_fraction' does not contain float datatype." in \
            err.value.stdout.decode('utf-8')


if __name__ == "__main__":
    exit(pytest.main([__file__]))
