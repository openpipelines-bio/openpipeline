import sys
import os
import h5py
import numpy as np
import pytest
import anndata as ad
import mudata as mu
import pandas as pd
import scipy.sparse as sp
from openpipeline_testutils.asserters import assert_annotation_objects_equal


## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

input_file = (
    f"{meta['resources_dir']}/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
)

raw_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_raw_feature_bc_matrix.h5"


@pytest.fixture(scope="module")
def empty_droplets(n_droplets=2000, min_counts=10, max_counts=100):
    # Real empty droplets from the raw (unfiltered) 10x h5 of the same run as the
    # input: low-count gene expression barcodes that are not called as cells.
    # Module-scoped, so the large raw file is only read once for all tests.
    with h5py.File(raw_file, "r") as f:
        matrix = f["matrix"]
        counts = sp.csc_matrix(
            (matrix["data"][:], matrix["indices"][:], matrix["indptr"][:]),
            shape=tuple(matrix["shape"][:]),
        )
        barcodes = matrix["barcodes"][:].astype(str)
        features = pd.DataFrame(
            {
                "gene_ids": matrix["features/id"][:].astype(str),
                "gene_symbol": matrix["features/name"][:].astype(str),
                "feature_types": matrix["features/feature_type"][:].astype(str),
            }
        ).set_index("gene_ids")

    is_gex = (features["feature_types"] == "Gene Expression").to_numpy()
    counts = counts[is_gex, :]
    library_size = np.asarray(counts.sum(axis=0)).ravel()

    cell_barcodes = mu.read_h5mu(input_file).mod["rna"].obs_names
    candidates = np.flatnonzero(
        (library_size >= min_counts)
        & (library_size <= max_counts)
        & ~np.isin(barcodes, cell_barcodes)
    )
    rng = np.random.default_rng(0)
    selected = np.sort(
        rng.choice(candidates, size=min(n_droplets, candidates.size), replace=False)
    )

    rna = ad.AnnData(
        X=counts[:, selected].T.tocsr(),
        obs=pd.DataFrame(index=barcodes[selected]),
        var=features[is_gex],
    )
    return mu.MuData({"rna": rna})


@pytest.fixture
def background_path(write_mudata_to_file, empty_droplets):
    return write_mudata_to_file(empty_droplets)


def assert_decontx_output(
    input_path,
    output_file,
    output_layer="decontx_counts",
    output_obs_contamination="decontx_contamination",
    output_obs_clusters="decontx_clusters",
):
    assert os.path.exists(output_file), "No output was created."

    output_mudata = mu.read_h5mu(output_file)
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

    # Delete output slots in order to test for equality with input data
    del output_rna.layers[output_layer]
    output_obs_columns = [output_obs_contamination, output_obs_clusters]
    output_rna.obs = output_rna.obs.drop(columns=output_obs_columns)
    output_mudata.update()
    assert_annotation_objects_equal(input_path, output_mudata)


def test_execution_without_background(run_component, random_h5mu_path):
    output_file = random_h5mu_path()

    logger.info("> Check whether decontx works without a background dataset")
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

    assert_decontx_output(input_file, output_file)


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
def background_with_batch_path(write_mudata_to_file, empty_droplets):
    # Copy so the shared module-scoped object is not modified
    bdata = empty_droplets.copy()
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


def test_execution_with_background(run_component, random_h5mu_path, background_path):
    output_file = random_h5mu_path()

    logger.info("> Check whether decontx works with a background dataset")
    run_component(
        [
            "--input",
            input_file,
            "--background",
            background_path,
            "--output",
            output_file,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(input_file, output_file)


def test_execution_with_cluster_label(
    run_component, random_h5mu_path, input_with_clusters_path
):
    output_file = random_h5mu_path()

    logger.info("> Check whether decontx works with cluster labels provided")
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

    assert_decontx_output(input_with_clusters_path, output_file)

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

    logger.info("> Check whether decontx works with batch labels provided")
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

    assert_decontx_output(input_with_batch_path, output_file)


def test_execution_with_batch_with_background(
    run_component, random_h5mu_path, input_with_batch_path, background_with_batch_path
):
    output_file = random_h5mu_path()

    logger.info(
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

    assert_decontx_output(input_with_batch_path, output_file)


def test_execution_without_delta_estimation(run_component, random_h5mu_path):
    output_default = random_h5mu_path()

    logger.info(
        "> Check whether decontx works with delta estimation enabled (baseline)"
    )
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

    logger.info("> Check whether decontx works with delta estimation disabled")
    run_component(
        [
            "--input",
            input_file,
            "--estimate_delta",
            "false",
            "--output",
            output_fixed_delta,
            "--output_compression",
            "gzip",
        ]
    )

    assert_decontx_output(input_file, output_default)
    assert_decontx_output(input_file, output_fixed_delta)

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

    logger.info("> Check whether decontx honors custom output layer/obs names")
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
        input_file,
        output_file,
        output_layer="custom_counts",
        output_obs_contamination="custom_contamination",
        output_obs_clusters="custom_clusters",
    )


def test_execution_with_input_layer(
    run_component, random_h5mu_path, input_with_layer_path
):
    output_default = random_h5mu_path()

    logger.info("> Check whether decontx works without --input_layer (baseline)")
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

    logger.info("> Check whether decontx works with --input_layer provided")
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

    assert_decontx_output(input_with_layer_path, output_default)
    assert_decontx_output(input_with_layer_path, output_layer)

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
