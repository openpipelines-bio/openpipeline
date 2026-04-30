import sys
import mudata as mu
import numpy as np
import pytest
import os

##VIASH START
par = {
    "input": "output.h5mu",
    "orig_input": [
        "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    ],
}

meta = {"resources_dir": "src/base"}
##VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

MIN_PERCENTILE = 0.05
MAX_PERCENTILE = 0.95
EXPECTED_RETENTION = MAX_PERCENTILE - MIN_PERCENTILE
RETENTION_TOLERANCE = 0.05


@pytest.fixture
def input_h5mu():
    return mu.read_h5mu(par["orig_input"][0])


@pytest.fixture
def output_h5mu():
    return mu.read_h5mu(par["input"])


@pytest.mark.parametrize("modality", ["rna", "prot"])
def test_filter_with_percentile_column_present(output_h5mu, modality):
    assert "filter_with_percentile" in output_h5mu[modality].obs.columns


@pytest.mark.parametrize("modality", ["rna", "prot"])
def test_total_counts_column_present(output_h5mu, modality):
    assert "total_counts" in output_h5mu[modality].obs.columns


@pytest.mark.parametrize("modality", ["rna", "prot"])
def test_cells_were_filtered(output_h5mu, input_h5mu, modality):
    assert output_h5mu[modality].n_obs < input_h5mu[modality].n_obs, (
        f"Expected percentile filter to remove cells from {modality}, "
        f"but output has {output_h5mu[modality].n_obs} cells vs input "
        f"{input_h5mu[modality].n_obs}."
    )


@pytest.mark.parametrize("modality", ["rna", "prot"])
def test_retention_matches_percentile_bounds(output_h5mu, input_h5mu, modality):
    retention = output_h5mu[modality].n_obs / input_h5mu[modality].n_obs
    assert abs(retention - EXPECTED_RETENTION) <= RETENTION_TOLERANCE, (
        f"Expected ~{EXPECTED_RETENTION:.0%} retention for {modality} "
        f"with percentiles [{MIN_PERCENTILE}, {MAX_PERCENTILE}], got {retention:.2%}."
    )


@pytest.mark.parametrize("modality", ["rna", "prot"])
def test_retained_cells_within_percentile_bounds(output_h5mu, input_h5mu, modality):
    # The filter is applied on log1p(total_counts) by default
    # (rna_log_transform_total_counts / prot_log_transform_total_counts default to true).
    input_total_counts = np.asarray(input_h5mu[modality].X.sum(axis=1)).flatten()
    log_input = np.log1p(input_total_counts)
    low, high = np.quantile(log_input, [MIN_PERCENTILE, MAX_PERCENTILE])

    output_total_counts = output_h5mu[modality].obs["total_counts"].to_numpy()
    log_output = np.log1p(output_total_counts)

    assert np.all((log_output >= low) & (log_output <= high)), (
        f"Retained {modality} cells have log1p(total_counts) outside "
        f"[{low:.3f}, {high:.3f}]: min={log_output.min():.3f}, "
        f"max={log_output.max():.3f}."
    )


if __name__ == "__main__":
    sys.exit(
        pytest.main(
            [
                "--import-mode=importlib",
                "-o",
                "python_files=script*.py .viash_script.py",
                os.path.dirname(__file__),
            ]
        )
    )
