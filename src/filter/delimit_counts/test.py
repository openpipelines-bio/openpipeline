import mudata as mu
import sys
import subprocess
import pytest
import os
import re

## VIASH START
meta = {
    "executable": "./target/executable/filter/delimit_counts/delimit_counts",
    "resources_dir": "resources_test/",
    "config": "/home/di/code/openpipeline/src/filter/delimit_counts/config.vsh.yaml",
}

## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def input_path():
    """Path to the pseudobulk test data"""
    return f"{meta['resources_dir']}/TS_Blood_filtered_pseudobulk.h5mu"


def test_simple_execution(run_component, tmp_path, input_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_path,
            "--obs_count_column",
            "n_cells",
            "--obs_name_filter",
            "filter_with_counts",
            "--min_count",
            "15",
            "--max_counts",
            "20",
            "--output",
            output_path,
        ]
    )

    assert os.path.exists(output_path), "Output file does not exist"
    adata_out = mu.read_h5ad(output_path, mod="rna")
    adata_in = mu.read_h5ad(input_path, mod="rna")

    adata_in.shape[0] > adata_out.shape[0], "Expected some cells to be filtered"
    adata_out.shape[1] == adata_in.shape[1], "Number of genes should be unchanged"
    adata_out.shape[0] == 2, "Expected 2 cells to remain after filtering"


def test_no_filtering(run_component, tmp_path, input_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_path,
            "--obs_count_column",
            "n_cells",
            "--obs_name_filter",
            "filter_with_counts",
            "--output",
            output_path,
        ]
    )

    assert os.path.exists(output_path), "Output file does not exist"
    adata_out = mu.read_h5ad(output_path, mod="rna")
    adata_in = mu.read_h5ad(input_path, mod="rna")

    adata_in.shape[0] == adata_out.shape[0], "Expected some cells to be filtered"
    adata_out.shape[1] == adata_in.shape[1], "Number of genes should be unchanged"
    adata_out.shape[0] == 16, "No filtering should have taken place"


def test_wrong_dtype(run_component, tmp_path, input_path):
    tmp_path = tmp_path / "wrong_dtype.h5mu"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--obs_count_column",
                "donor_id",
                "--obs_name_filter",
                "filter_with_counts",
                "--output",
                tmp_path,
            ]
        )

    assert re.search(
        r"Column 'donor_id' does not contain integer datatype.",
        err.value.stdout.decode("utf-8"),
    ), f"Expected error message not found: {err.value.stdout.decode('utf-8')}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
