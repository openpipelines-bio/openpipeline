import mudata as mu
import sys
from pathlib import Path
import pytest
import numpy as np
import pandas as pd
import subprocess
import re

## VIASH START
meta = {
    "executable": "./target/executable/filter/filter_with_quantile/filter_with_quantile",
    "resources_dir": "src/utils",
    "config": "./src/filter/filter_with_quantile/config.vsh.yaml",
}

## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def base_h5mu():
    """Load the base h5mu data without additional columns."""
    path = f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
    return mu.read_h5mu(path)


@pytest.fixture
def input_h5mu(base_h5mu):
    """Base h5mu data with random normal columns for standard tests."""
    input_data = base_h5mu.copy()

    # Generate normally distributed integers between 0 and 100 for .obs
    random_normal = np.random.normal(50, 5, input_data.mod["rna"].n_obs)
    input_data.mod["rna"].obs["random_normal_int"] = np.clip(
        random_normal, 0, 100
    ).astype(int)

    # Generate normally distributed integers between 0 and 100 for .var
    random_normal_var = np.random.normal(50, 5, input_data.mod["rna"].n_vars)
    input_data.mod["rna"].var["random_normal_int"] = np.clip(
        random_normal_var, 0, 100
    ).astype(int)

    return input_data


@pytest.fixture
def input_h5mu_with_nan(input_h5mu):
    """Factory fixture for H5mu data with different NaN types for error testing."""

    def _create_data_with_nan(na_value):
        input_data = input_h5mu.copy()

        # Add columns with NaN values for testing error handling
        # Use appropriate dtype based on NA type
        obs_values = np.random.normal(50, 5, input_data.mod["rna"].n_obs)
        var_values = np.random.normal(50, 5, input_data.mod["rna"].n_vars)

        # Add NA values to numerical data
        obs_series = pd.Series(obs_values, dtype="float64")
        var_series = pd.Series(var_values, dtype="float64")

        obs_series.iloc[0:3] = na_value  # Add specific NA values
        var_series.iloc[0:2] = na_value  # Add specific NA values

        input_data.mod["rna"].obs["random_with_nan"] = obs_series
        input_data.mod["rna"].var["random_with_nan"] = var_series

        return input_data

    return _create_data_with_nan


@pytest.fixture
def input_h5mu_skewed(base_h5mu):
    """H5mu data with skewed distribution that benefits from log1p transformation."""
    input_data = base_h5mu.copy()

    # Generate exponentially distributed (skewed) data that will benefit from log1p transformation
    # Using exponential distribution with different scales for obs and var
    exp_obs_values = np.random.exponential(scale=10, size=input_data.mod["rna"].n_obs)
    exp_var_values = np.random.exponential(scale=5, size=input_data.mod["rna"].n_vars)

    # Convert to integers to simulate count-like data
    input_data.mod["rna"].obs["skewed_counts"] = exp_obs_values.astype(int)
    input_data.mod["rna"].var["skewed_counts"] = exp_var_values.astype(int)

    return input_data


@pytest.fixture
def input_path_skewed(input_h5mu_skewed, random_h5mu_path):
    path = random_h5mu_path()
    input_h5mu_skewed.write(path)
    return path


@pytest.fixture
def input_path(input_h5mu, random_h5mu_path):
    path = random_h5mu_path()
    input_h5mu.write(path)
    return path


@pytest.fixture
def input_path_with_nan_factory(input_h5mu_with_nan, random_h5mu_path):
    """Factory fixture for creating input paths with different NaN types."""

    def _create_path_with_nan(na_value):
        data = input_h5mu_with_nan(na_value)
        path = random_h5mu_path()
        data.write(path)
        return path

    return _create_path_with_nan


@pytest.fixture
def input_n_rna_obs(input_h5mu):
    return input_h5mu.mod["rna"].n_obs


@pytest.fixture
def input_n_prot_obs(input_h5mu):
    return input_h5mu.mod["prot"].n_obs


@pytest.fixture
def input_n_rna_vars(input_h5mu):
    return input_h5mu.mod["rna"].n_vars


@pytest.fixture
def input_n_prot_vars(input_h5mu):
    return input_h5mu.mod["prot"].n_vars


def test_filter_mask(
    run_component,
    input_path,
    random_h5mu_path,
    input_n_rna_obs,
    input_n_prot_obs,
    input_n_rna_vars,
    input_n_prot_vars,
):
    output_path = random_h5mu_path()

    args = [
        "--input",
        input_path,
        "--output",
        output_path,
        "--obs_column",
        "random_normal_int",
        "--var_column",
        "random_normal_int",
        "--obs_min_quantile",
        "0.1",
        "--obs_max_quantile",
        "0.9",
        "--var_min_quantile",
        "0.05",
        "--var_max_quantile",
        "0.95",
    ]

    run_component(args)

    assert Path(output_path).is_file()
    mu_out = mu.read_h5mu(output_path)
    rna_mod = mu_out.mod["rna"]
    prot_mod = mu_out.mod["prot"]

    assert "filter_with_quantile" in rna_mod.obs
    assert "filter_with_quantile" in rna_mod.var

    new_obs = rna_mod.n_obs
    new_vars = rna_mod.n_vars
    assert new_obs == input_n_rna_obs
    assert new_vars == input_n_rna_vars

    assert prot_mod.n_obs == input_n_prot_obs
    assert prot_mod.n_vars == input_n_prot_vars

    obs_filter = rna_mod.obs["filter_with_quantile"]
    obs_filter_fraction = sum(obs_filter) / len(obs_filter)
    var_filter = rna_mod.var["filter_with_quantile"]
    var_filter_fraction = sum(var_filter) / len(var_filter)

    assert np.isclose(obs_filter_fraction, 0.8, atol=0.1), (
        f"Expected ~80% of obs to be kept, but got {obs_filter_fraction:.2f}"
    )
    assert np.isclose(var_filter_fraction, 0.9, atol=0.05), (
        f"Expected ~90% of vars to be kept, but got {var_filter_fraction:.2f}"
    )


def test_subset(
    run_component,
    input_path,
    random_h5mu_path,
    input_n_rna_obs,
    input_n_prot_obs,
    input_n_rna_vars,
    input_n_prot_vars,
):
    output_path = random_h5mu_path()

    args = [
        "--input",
        input_path,
        "--output",
        output_path,
        "--obs_column",
        "random_normal_int",
        "--var_column",
        "random_normal_int",
        "--obs_min_quantile",
        "0.1",
        "--obs_max_quantile",
        "0.9",
        "--var_min_quantile",
        "0.05",
        "--var_max_quantile",
        "0.95",
        "--do_subset",
        "True",
    ]

    run_component(args)

    assert Path(output_path).is_file()
    mu_out = mu.read_h5mu(output_path)
    rna_mod = mu_out.mod["rna"]
    prot_mod = mu_out.mod["prot"]

    assert "filter_with_quantile" in rna_mod.obs
    assert "filter_with_quantile" in rna_mod.var

    new_obs = rna_mod.n_obs
    new_vars = rna_mod.n_vars

    assert np.isclose(new_obs / input_n_rna_obs, 0.8, atol=0.1), (
        "Expected .obs to be subset."
    )
    assert np.isclose(new_vars / input_n_rna_vars, 0.9, atol=0.05), (
        "Expected .var to be subset."
    )

    assert prot_mod.n_obs == input_n_prot_obs
    assert prot_mod.n_vars == input_n_prot_vars


def test_no_quantiles(
    run_component, input_path, random_h5mu_path, input_n_rna_obs, input_n_rna_vars
):
    output_path = random_h5mu_path()

    args = [
        "--input",
        input_path,
        "--output",
        output_path,
        "--obs_column",
        "random_normal_int",
        "--var_column",
        "random_normal_int",
    ]

    run_component(args)

    assert Path(output_path).is_file()
    mu_out = mu.read_h5mu(output_path)
    rna_mod = mu_out.mod["rna"]

    assert "filter_with_quantile" not in rna_mod.obs
    assert "filter_with_quantile" not in rna_mod.var


def test_no_filter_columns(
    run_component, input_path, random_h5mu_path, input_n_rna_obs, input_n_rna_vars
):
    output_path = random_h5mu_path()

    args = ["--input", input_path, "--output", output_path]

    run_component(args)

    assert Path(output_path).is_file()
    mu_out = mu.read_h5mu(output_path)
    rna_mod = mu_out.mod["rna"]

    assert "filter_with_quantile" not in rna_mod.obs
    assert "filter_with_quantile" not in rna_mod.var

    assert rna_mod.n_obs == input_n_rna_obs, (
        "Expected no filtering to be applied when no filter columns are set."
    )
    assert rna_mod.n_vars == input_n_rna_vars, (
        "Expected no filtering to be applied when no filter columns are set."
    )


def test_raises_with_non_numeric_column(
    run_component,
    input_path,
    random_h5mu_path,
):
    output_path = random_h5mu_path()
    # fails because input data are not correctly lognormalized
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                output_path,
                "--var_column",
                "gene_symbol",
                "--var_min_quantile",
                "0.1",
            ]
        )
    assert re.search(
        r"Column 'gene_symbol' must contain numeric data for quantile filtering",
        err.value.stdout.decode("utf-8"),
    )


def test_raises_with_non_existent_column(
    run_component,
    input_path,
    random_h5mu_path,
):
    output_path = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                output_path,
                "--var_column",
                "non_existent_column",
            ]
        )
    assert re.search(
        r"Column 'non_existent_column' not found in .var. Available columns: .*",
        err.value.stdout.decode("utf-8"),
    )

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                output_path,
                "--obs_column",
                "non_existent_column",
            ]
        )
    assert re.search(
        r"Column 'non_existent_column' not found in .obs. Available columns: .*",
        err.value.stdout.decode("utf-8"),
    )


@pytest.mark.parametrize("na_value,na_name", [(np.nan, "np.nan"), (pd.NA, "pd.NA")])
def test_raises_with_nan(
    run_component,
    input_path_with_nan_factory,
    random_h5mu_path,
    na_value,
    na_name,
):
    # Create input path with specific NA value
    input_path = input_path_with_nan_factory(na_value)
    output_path = random_h5mu_path()

    # Test with obs column containing NaN values
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                output_path,
                "--obs_column",
                "random_with_nan",
                "--obs_min_quantile",
                "0.1",
            ]
        )

    error_message = err.value.stdout.decode("utf-8")
    assert re.search(
        r"Column contains NaN values. Please clean the data before applying quantile filtering.",
        error_message,
    ), (
        f"Expected NaN error message for {na_name} in obs column, but got: {error_message}"
    )

    # Test with var column containing NaN values
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                output_path,
                "--var_column",
                "random_with_nan",
                "--var_min_quantile",
                "0.1",
            ]
        )

    error_message = err.value.stdout.decode("utf-8")
    assert re.search(
        r"Column contains NaN values. Please clean the data before applying quantile filtering.",
        error_message,
    ), (
        f"Expected NaN error message for {na_name} in var column, but got: {error_message}"
    )


def test_log1p_transformation(
    run_component,
    input_path_skewed,
    random_h5mu_path,
    input_h5mu_skewed,
):
    """Test that log1p transformation is applied correctly and creates expected columns."""
    output_path = random_h5mu_path()

    args = [
        "--input",
        input_path_skewed,
        "--output",
        output_path,
        "--obs_column",
        "skewed_counts",
        "--var_column",
        "skewed_counts",
        "--obs_log1p_transform",
        "True",
        "--var_log1p_transform",
        "True",
        "--obs_min_quantile",
        "0.1",
        "--obs_max_quantile",
        "0.9",
        "--var_min_quantile",
        "0.05",
        "--var_max_quantile",
        "0.95",
    ]

    run_component(args)

    assert Path(output_path).is_file()
    mu_out = mu.read_h5mu(output_path)
    rna_mod = mu_out.mod["rna"]

    # Check that log1p transformed columns were created
    assert "log1p_skewed_counts" in rna_mod.obs.columns, (
        "Expected log1p_skewed_counts column to be created in .obs"
    )
    assert "log1p_skewed_counts" in rna_mod.var.columns, (
        "Expected log1p_skewed_counts column to be created in .var"
    )

    # Check that filter columns were created
    assert "filter_with_quantile" in rna_mod.obs
    assert "filter_with_quantile" in rna_mod.var

    # Validate log1p transformation correctness
    original_obs_values = input_h5mu_skewed.mod["rna"].obs["skewed_counts"].values
    original_var_values = input_h5mu_skewed.mod["rna"].var["skewed_counts"].values

    expected_obs_log1p = np.log1p(original_obs_values)
    expected_var_log1p = np.log1p(original_var_values)

    actual_obs_log1p = rna_mod.obs["log1p_skewed_counts"].values
    actual_var_log1p = rna_mod.var["log1p_skewed_counts"].values

    # Check that log1p transformation was applied correctly
    np.testing.assert_array_almost_equal(
        expected_obs_log1p,
        actual_obs_log1p,
        decimal=10,
        err_msg="log1p transformation not applied correctly to obs column",
    )
    np.testing.assert_array_almost_equal(
        expected_var_log1p,
        actual_var_log1p,
        decimal=10,
        err_msg="log1p transformation not applied correctly to var column",
    )

    # Check that filtering was applied (some values should be filtered out)
    obs_filter = rna_mod.obs["filter_with_quantile"]
    var_filter = rna_mod.var["filter_with_quantile"]

    obs_filter_fraction = sum(obs_filter) / len(obs_filter)
    var_filter_fraction = sum(var_filter) / len(var_filter)

    # Should be approximately 80% and 90% respectively, allowing some tolerance
    assert np.isclose(obs_filter_fraction, 0.8, atol=0.5), (
        f"Expected ~80% of obs to be kept after log1p transformation and filtering, "
        f"but got {obs_filter_fraction:.2f}"
    )
    assert np.isclose(var_filter_fraction, 0.9, atol=0.1), (
        f"Expected ~90% of vars to be kept after log1p transformation and filtering, "
        f"but got {var_filter_fraction:.2f}"
    )


def test_log1p_custom_column_names(
    run_component,
    input_path_skewed,
    random_h5mu_path,
    input_h5mu_skewed,
):
    """Test that custom log1p column names work correctly."""
    output_path = random_h5mu_path()

    custom_obs_column = "custom_obs_log_transformed"
    custom_var_column = "custom_var_log_transformed"

    args = [
        "--input",
        input_path_skewed,
        "--output",
        output_path,
        "--obs_column",
        "skewed_counts",
        "--var_column",
        "skewed_counts",
        "--obs_log1p_transform",
        "True",
        "--var_log1p_transform",
        "True",
        "--obs_log1p_column",
        custom_obs_column,
        "--var_log1p_column",
        custom_var_column,
        "--obs_min_quantile",
        "0.1",
        "--obs_max_quantile",
        "0.9",
        "--var_min_quantile",
        "0.05",
        "--var_max_quantile",
        "0.95",
    ]

    run_component(args)

    assert Path(output_path).is_file()
    mu_out = mu.read_h5mu(output_path)
    rna_mod = mu_out.mod["rna"]

    # Check that custom named log1p columns were created
    assert custom_obs_column in rna_mod.obs.columns, (
        f"Expected {custom_obs_column} column to be created in .obs"
    )
    assert custom_var_column in rna_mod.var.columns, (
        f"Expected {custom_var_column} column to be created in .var"
    )

    # Make sure default column names were NOT created
    assert "log1p_skewed_counts" not in rna_mod.obs.columns, (
        "Default log1p_skewed_counts should not be created when custom name is provided"
    )
    assert "log1p_skewed_counts" not in rna_mod.var.columns, (
        "Default log1p_skewed_counts should not be created when custom name is provided"
    )

    # Validate transformation correctness with custom column names
    original_obs_values = input_h5mu_skewed.mod["rna"].obs["skewed_counts"].values
    original_var_values = input_h5mu_skewed.mod["rna"].var["skewed_counts"].values

    expected_obs_log1p = np.log1p(original_obs_values)
    expected_var_log1p = np.log1p(original_var_values)

    actual_obs_log1p = rna_mod.obs[custom_obs_column].values
    actual_var_log1p = rna_mod.var[custom_var_column].values

    # Check that log1p transformation was applied correctly to custom columns
    np.testing.assert_array_almost_equal(
        expected_obs_log1p,
        actual_obs_log1p,
        decimal=10,
        err_msg="log1p transformation not applied correctly to custom obs column",
    )
    np.testing.assert_array_almost_equal(
        expected_var_log1p,
        actual_var_log1p,
        decimal=10,
        err_msg="log1p transformation not applied correctly to custom var column",
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
