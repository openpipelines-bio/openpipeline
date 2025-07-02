import os
import subprocess
import mudata as mu
import sys
import pytest
import re
import numpy as np
import pandas as pd


## VIASH START
meta = {
    "resources_dir": "resources_test/"
}
## VIASH END

sys.path.append(meta["resources_dir"])


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/TS_Blood_filtered.h5mu"


@pytest.fixture
def input_data(input_path):
    mu_in = mu.read_h5mu(input_path)
    return mu_in


@pytest.fixture
def annotated_test_data(input_data):
    rna_in = input_data.mod["rna"]
    np.random.seed(0)
    n_cells = rna_in.n_obs
    treatment = np.random.choice(['ctrl', 'stim'], size=n_cells, p=[0.5, 0.5])
    disease = np.random.choice(['healthy', 'diseased'], size=n_cells, p=[0.5, 0.5])
    rna_in.obs['treatment'] = treatment
    rna_in.obs['disease'] = disease
    return input_data


@pytest.fixture
def annotated_test_data_path(tmp_path, annotated_test_data):
    temp_h5mu = tmp_path / "de_test_data.h5mu"
    annotated_test_data.write_h5mu(temp_h5mu)
    return temp_h5mu


def test_simple_execution(run_component, random_h5mu_path, annotated_test_data_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            annotated_test_data_path,
            "--output",
            output_path,
            "--obs_grouping",
            "cell_type",
            "--obs_sample_conditions",
            "treatment",
            "--output_compression",
            "gzip",
        ]
    )

    assert os.path.exists(output_path), "Output query file does not exist"
    mdata = mu.read_h5mu(output_path)
    adata = mdata.mod["rna"]

    expected_obs = ["treatment", "cell_type", "n_cells"]
    assert all(col in adata.obs for col in expected_obs), f"Expected columns {expected_obs} not found in .obs"
    assert adata.shape[0] == 8, "Expected a total of 8 pseudobulk samples in the output"


def test_multiple_factors(run_component, random_h5mu_path, annotated_test_data_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            annotated_test_data_path,
            "--output",
            output_path,
            "--obs_grouping",
            "cell_type",
            "--obs_sample_conditions",
            "disease",
            "--obs_sample_conditions",
            "donor_id",
            "--obs_sample_conditions",
            "treatment",
            "--min_num_cells_per_sample",
            "5",
            "--output_compression",
            "gzip",
        ]
    )

    assert os.path.exists(output_path), "Output query file does not exist"
    mdata = mu.read_h5mu(output_path)
    adata = mdata.mod["rna"]

    expected_obs = ["treatment", "cell_type", "n_cells"]
    assert all(col in adata.obs for col in expected_obs), f"Expected columns {expected_obs} not found in .obs"
    assert adata.shape[0] == 16, "Expected a total of 16 pseudobulk samples in the output"


def test_pseudo_replicates(run_component, random_h5mu_path, annotated_test_data_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            annotated_test_data_path,
            "--output",
            output_path,
            "--obs_grouping",
            "cell_type",
            "--obs_sample_conditions",
            "treatment",
            "--pseudo_replicates",
            "2",
            "--output_compression",
            "gzip",
        ]
    )

    assert os.path.exists(output_path), "Output query file does not exist"
    mdata = mu.read_h5mu(output_path)
    adata = mdata.mod["rna"]

    expected_obs = ["treatment", "cell_type", "n_cells"]
    assert all(col in adata.obs for col in expected_obs), f"Expected columns {expected_obs} not found in .obs"
    assert adata.shape[0] == 16, "Expected a total of 8 pseudobulk samples in the output"


def test_filtering(run_component, random_h5mu_path, annotated_test_data_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            annotated_test_data_path,
            "--output",
            output_path,
            "--obs_grouping",
            "cell_type",
            "--obs_sample_conditions",
            "treatment",
            "--min_num_cells_per_sample",
            "50",
            "--output_compression",
            "gzip",
        ]
    )

    assert os.path.exists(output_path), "Output query file does not exist"
    mdata = mu.read_h5mu(output_path)
    adata = mdata.mod["rna"]

    expected_obs = ["treatment", "cell_type", "n_cells"]
    assert all(col in adata.obs for col in expected_obs), f"Expected columns {expected_obs} not found in .obs"
    assert adata.shape[0] == 4, "Expected a total of 8 pseudobulk samples in the output"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
