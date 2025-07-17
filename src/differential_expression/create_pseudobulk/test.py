import os
import mudata as mu
import sys
import pytest

## VIASH START
meta = {"resources_dir": "resources_test/"}
## VIASH END

sys.path.append(meta["resources_dir"])

input_path = meta["resources_dir"] + "/TS_Blood_filtered.h5mu"


def test_simple_execution(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
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
    assert all(col in adata.obs for col in expected_obs), (
        f"Expected columns {expected_obs} not found in .obs"
    )
    assert adata.shape[0] == 8, "Expected a total of 8 pseudobulk samples in the output"


def test_multiple_factors(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
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
    assert all(col in adata.obs for col in expected_obs), (
        f"Expected columns {expected_obs} not found in .obs"
    )
    assert adata.shape[0] == 16, (
        "Expected a total of 16 pseudobulk samples in the output"
    )


def test_pseudo_replicates(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
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
    assert all(col in adata.obs for col in expected_obs), (
        f"Expected columns {expected_obs} not found in .obs"
    )
    assert adata.shape[0] == 16, (
        "Expected a total of 8 pseudobulk samples in the output"
    )


def test_filtering(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
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
    assert all(col in adata.obs for col in expected_obs), (
        f"Expected columns {expected_obs} not found in .obs"
    )
    assert adata.shape[0] == 4, "Expected a total of 8 pseudobulk samples in the output"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
