import mudata as mu
import numpy as np
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


@pytest.fixture
def var_input_path(tmp_path, input_path):
    """Test data with an added numeric .var column to threshold on."""
    output_path = tmp_path / "with_var_column.h5mu"
    mdata = mu.read_h5mu(input_path)
    adata = mdata.mod["rna"]
    adata.var["var_counts"] = np.arange(adata.n_vars) % 10
    mdata.write(output_path)
    return output_path


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
            "--max_count",
            "20",
            "--do_subset",
            "True",
            "--output",
            output_path,
        ]
    )

    assert os.path.exists(output_path), "Output file does not exist"
    adata_out = mu.read_h5ad(output_path, mod="rna")
    adata_in = mu.read_h5ad(input_path, mod="rna")

    assert adata_out.shape[0] == 2, "Expected 2 cells to remain after filtering"
    assert adata_out.shape[1] == adata_in.shape[1], (
        "Number of genes should be unchanged"
    )


def test_no_filtering(run_component, tmp_path, input_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            input_path,
            "--obs_count_column",
            "n_cells",
            "--do_subset",
            "True",
            "--obs_name_filter",
            "filter_with_counts",
            "--output",
            output_path,
        ]
    )

    assert os.path.exists(output_path), "Output file does not exist"
    adata_out = mu.read_h5ad(output_path, mod="rna")
    adata_in = mu.read_h5ad(input_path, mod="rna")

    assert adata_out.shape[0] == adata_in.shape[0], (
        "No filtering should have taken place"
    )
    assert adata_out.shape[1] == adata_in.shape[1], (
        "Number of genes should be unchanged"
    )


def test_no_subset(run_component, tmp_path, input_path):
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
            "--max_count",
            "20",
            "--output",
            output_path,
        ]
    )

    assert os.path.exists(output_path), "Output file does not exist"
    adata_out = mu.read_h5ad(output_path, mod="rna")
    adata_in = mu.read_h5ad(input_path, mod="rna")

    assert adata_out.shape[0] == adata_in.shape[0], "No cells should be removed"
    assert "filter_with_counts" in adata_out.obs, "Filter column should be stored"
    assert adata_out.obs["filter_with_counts"].sum() == 2, (
        "Expected 2 cells to be retained by the filter mask"
    )


def test_var_filtering(run_component, tmp_path, var_input_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            var_input_path,
            "--var_count_column",
            "var_counts",
            "--var_name_filter",
            "filter_var_counts",
            "--min_count",
            "5",
            "--do_subset",
            "True",
            "--output",
            output_path,
        ]
    )

    assert os.path.exists(output_path), "Output file does not exist"
    adata_out = mu.read_h5ad(output_path, mod="rna")
    adata_in = mu.read_h5ad(var_input_path, mod="rna")

    expected_vars = int((adata_in.var["var_counts"] >= 5).sum())
    assert adata_out.shape[1] == expected_vars, (
        "Genes below the threshold should be removed"
    )
    assert adata_out.shape[0] == adata_in.shape[0], (
        "Number of cells should be unchanged"
    )


def test_obs_and_var_filtering(run_component, tmp_path, var_input_path):
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            var_input_path,
            "--obs_count_column",
            "n_cells",
            "--obs_name_filter",
            "filter_with_counts",
            "--var_count_column",
            "var_counts",
            "--var_name_filter",
            "filter_var_counts",
            "--min_count",
            "5",
            "--output",
            output_path,
        ]
    )

    assert os.path.exists(output_path), "Output file does not exist"
    adata_out = mu.read_h5ad(output_path, mod="rna")
    assert "filter_with_counts" in adata_out.obs, "obs filter column should be stored"
    assert "filter_var_counts" in adata_out.var, "var filter column should be stored"


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
        r"Column 'donor_id' does not contain numeric datatype.",
        err.value.stdout.decode("utf-8"),
    ), f"Expected error message not found: {err.value.stdout.decode('utf-8')}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
