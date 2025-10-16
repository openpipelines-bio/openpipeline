import os
import mudata as mu
import sys
import re
import subprocess
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
            "--obs_label",
            "cell_type",
            "--obs_groups",
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


def test_log_normalized_counts(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                output_path,
                "--input_layer",
                "log_normalized",
                "--obs_label",
                "cell_type",
                "--obs_groups",
                "treatment",
                "--output_compression",
                "gzip",
            ]
        )

    assert re.search(
        r"ValueError: Input layer must contain raw counts.",
        err.value.stdout.decode("utf-8"),
    )


def test_multiple_factors(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--obs_label",
            "cell_type",
            "--obs_groups",
            "disease",
            "--obs_groups",
            "donor_id",
            "--obs_groups",
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
    assert adata.shape[0] == 16, (
        "Expected a total of 16 pseudobulk samples in the output"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
