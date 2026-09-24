import sys
import os
import numpy as np
import scipy.sparse as sp
import pytest
import mudata as mu
import pandas as pd
from openpipeline_testutils.asserters import assert_annotation_objects_equal


def _assert_matrices_equal(a, b, msg):
    if sp.issparse(a) or sp.issparse(b):
        a, b = sp.csr_matrix(a), sp.csr_matrix(b)
        assert a.shape == b.shape and (a != b).nnz == 0, msg
    else:
        assert np.array_equal(np.asarray(a), np.asarray(b)), msg


## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

input_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)

background_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5mu"
)


def assert_decontx_output(
    output_file,
    output_layer="decontx_counts",
    output_obs_contamination="decontx_contamination",
    output_obs_clusters="decontx_clusters",
):
    assert os.path.exists(output_file), "No output was created."

    input_mudata = mu.read_h5mu(input_file)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    output_rna = output_mudata.mod["rna"]
    assert output_layer in output_rna.layers, (
        f"Output should contain a '{output_layer}' layer."
    )
    assert output_obs_contamination in output_rna.obs, (
        f"Output should contain a '{output_obs_contamination}' .obs column."
    )
    assert output_obs_clusters in output_rna.obs, (
        f"Output should contain a '{output_obs_clusters}' .obs column."
    )
    assert output_rna.obs[output_obs_contamination].between(0, 1).all(), (
        "Contamination values should be between 0 and 1."
    )

    input_rna = input_mudata.mod["rna"]
    _assert_matrices_equal(
        input_rna.X, output_rna.X, "rna.X should be unchanged by the component."
    )

    for layer_name, layer_matrix in input_rna.layers.items():
        assert layer_name in output_rna.layers, (
            f"Pre-existing layer '{layer_name}' should still be present in the output."
        )
        _assert_matrices_equal(
            layer_matrix,
            output_rna.layers[layer_name],
            f"Pre-existing layer '{layer_name}' should be unchanged by the component.",
        )


def test_execution_without_background(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    print("> Check whether decontx works without a background dataset")
    run_component(
        [
            "--input",
            input_file,
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_file)


@pytest.fixture
def input_with_clusters_path(write_mudata_to_file):
    mdata = mu.read_h5mu(input_file)
    rna = mdata.mod["rna"]
    rna.obs["cell_cluster"] = pd.Categorical(
        np.random.default_rng(0).integers(0, 3, size=rna.n_obs).astype(str)
    )
    return write_mudata_to_file(mdata)


@pytest.fixture
def input_with_batch_path(write_mudata_to_file):
    mdata = mu.read_h5mu(input_file)
    rna = mdata.mod["rna"]
    rna.obs["batch"] = pd.Categorical(
        np.random.default_rng(0).integers(0, 2, size=rna.n_obs).astype(str)
    )
    return write_mudata_to_file(mdata)


@pytest.fixture
def background_with_batch_path(write_mudata_to_file):
    bdata = mu.read_h5mu(background_file)
    rna = bdata.mod["rna"]
    rna.obs["batch"] = pd.Categorical(
        np.random.default_rng(0).integers(0, 2, size=rna.n_obs).astype(str)
    )
    return write_mudata_to_file(bdata)


@pytest.fixture
def input_with_layer_path(write_mudata_to_file):
    mdata = mu.read_h5mu(input_file)
    rna = mdata.mod["rna"]
    # Permute cells (not just scale) so the layer's per-cell profiles genuinely
    # differ from .X, rather than merely differing by a scale-invariant factor
    # that DecontX's proportion-based estimate could shrug off.
    perm = np.random.default_rng(0).permutation(rna.n_obs)
    rna.layers["shuffled_counts"] = rna.X[perm, :]
    return write_mudata_to_file(mdata)


def test_execution_with_background(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    print("> Check whether decontx works with a background dataset")
    run_component(
        [
            "--input",
            input_file,
            "--background",
            background_file,
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_file)


def test_execution_with_cluster_label(
    run_component, random_h5mu_path, input_with_clusters_path
):
    output_file = random_h5mu_path()

    print("> Check whether decontx works with cluster labels provided")
    run_component(
        [
            "--input",
            input_with_clusters_path,
            "--input_obs_clusters",
            "cell_cluster",
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_file)

    input_clusters = (
        mu.read_h5mu(input_with_clusters_path)
        .mod["rna"]
        .obs["cell_cluster"]
        .astype(str)
    )
    output_clusters = (
        mu.read_h5mu(output_file).mod["rna"].obs["decontx_clusters"].astype(str)
    )

    # decontX re-encodes cluster labels as a 1-indexed integer factor, so it
    # won't echo the original label spellings back - only the grouping of
    # cells into clusters is preserved (a bijection between input and output
    # labels), not the label values themselves.
    grouping = pd.DataFrame(
        {"input": input_clusters.values, "output": output_clusters.values}
    )
    assert grouping.groupby("input")["output"].nunique().eq(1).all(), (
        "Cells sharing an --input_obs_clusters label should map to a single decontx_clusters label."
    )
    assert grouping["output"].nunique() == grouping["input"].nunique(), (
        "decontx_clusters should use a distinct label for each --input_obs_clusters group."
    )


def test_execution_with_batch(run_component, random_h5mu_path, input_with_batch_path):
    output_file = random_h5mu_path()

    print("> Check whether decontx works with batch labels provided")
    run_component(
        [
            "--input",
            input_with_batch_path,
            "--input_obs_batch",
            "batch",
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_file)


def test_execution_with_batch_with_background(
    run_component, random_h5mu_path, input_with_batch_path, background_with_batch_path
):
    output_file = random_h5mu_path()

    print(
        "> Check whether decontx works with batch labels and a background dataset provided"
    )
    run_component(
        [
            "--input",
            input_with_batch_path,
            "--input_obs_batch",
            "batch",
            "--background",
            background_with_batch_path,
            "--background_obs_batch",
            "batch",
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_file)


def test_execution_without_delta_estimation(run_component, random_h5mu_path):
    output_default = random_h5mu_path()

    print("> Check whether decontx works with delta estimation enabled (baseline)")
    run_component(
        [
            "--input",
            input_file,
            "--output",
            output_default,
            "--output_compression",
            "gzip",
        ]
    )

    output_fixed_delta = random_h5mu_path()

    print("> Check whether decontx works with delta estimation disabled")
    run_component(
        [
            "--input",
            input_file,
            "--estimate_delta",
            "--output",
            output_fixed_delta,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_fixed_delta)

    contamination_default = (
        mu.read_h5mu(output_default).mod["rna"].obs["decontx_contamination"]
    )
    contamination_fixed_delta = (
        mu.read_h5mu(output_fixed_delta).mod["rna"].obs["decontx_contamination"]
    )
    assert not np.allclose(
        contamination_default.values, contamination_fixed_delta.values
    ), (
        "Disabling delta re-estimation should change DecontX's converged result; "
        "identical results suggest --estimate_delta was silently ignored."
    )


def test_execution_with_custom_output_names(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    print("> Check whether decontx honors custom output layer/obs names")
    run_component(
        [
            "--input",
            input_file,
            "--output",
            output_file,
            "--output_compression",
            "gzip",
            "--output_layer",
            "custom_counts",
            "--output_obs_contamination",
            "custom_contamination",
            "--output_obs_clusters",
            "custom_clusters",
        ]
    )

    assert_decontx_output(
        output_file,
        output_layer="custom_counts",
        output_obs_contamination="custom_contamination",
        output_obs_clusters="custom_clusters",
    )

    output_rna = mu.read_h5mu(output_file).mod["rna"]
    assert "decontx_counts" not in output_rna.layers, (
        "Default 'decontx_counts' layer should not be written when --output_layer is set."
    )
    assert "decontx_contamination" not in output_rna.obs, (
        "Default 'decontx_contamination' column should not be written when "
        "--output_obs_contamination is set."
    )
    assert "decontx_clusters" not in output_rna.obs, (
        "Default 'decontx_clusters' column should not be written when "
        "--output_obs_clusters is set."
    )


def test_execution_with_input_layer(
    run_component, random_h5mu_path, input_with_layer_path
):
    output_default = random_h5mu_path()

    print("> Check whether decontx works without --input_layer (baseline)")
    run_component(
        [
            "--input",
            input_with_layer_path,
            "--output",
            output_default,
            "--output_compression",
            "gzip",
        ]
    )

    output_layer = random_h5mu_path()

    print("> Check whether decontx works with --input_layer provided")
    run_component(
        [
            "--input",
            input_with_layer_path,
            "--input_layer",
            "shuffled_counts",
            "--output",
            output_layer,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(output_layer)

    input_layer_data = (
        mu.read_h5mu(input_with_layer_path).mod["rna"].layers["shuffled_counts"]
    )
    output_layer_data = mu.read_h5mu(output_layer).mod["rna"].layers["shuffled_counts"]
    _assert_matrices_equal(
        input_layer_data,
        output_layer_data,
        "Pre-existing 'shuffled_counts' layer should be unchanged by the component.",
    )

    contamination_default = (
        mu.read_h5mu(output_default).mod["rna"].obs["decontx_contamination"]
    )
    contamination_layer = (
        mu.read_h5mu(output_layer).mod["rna"].obs["decontx_contamination"]
    )
    assert not np.allclose(contamination_default.values, contamination_layer.values), (
        "Using --input_layer should change DecontX's input and therefore its output; "
        "identical results suggest the component silently fell back to .X."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
