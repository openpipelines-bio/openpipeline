import os
import subprocess
import scanpy as sc
import mudata as mu
import sys
import pytest
import re
import pandas as pd
import numpy as np
import anndata as ad

from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "resources_dir": "resources_test/",
    "config": "./src/feature_annotation/highly_variable_features_scanpy/config.vsh.yaml",
    "executable": "./target/executable/feature_annotation/highly_variable_features_scanpy/highly_variable_features_scanpy",
}
## VIASH END

sys.path.append(meta["resources_dir"])


@pytest.fixture
def input_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def input_data(input_path):
    mu_in = mu.read_h5mu(input_path)
    return mu_in


@pytest.fixture
def lognormed_test_data(input_data):
    rna_in = input_data.mod["rna"]
    assert "filter_with_hvg" not in rna_in.var.columns
    log_transformed = sc.pp.log1p(rna_in, copy=True)
    rna_in.layers["log_transformed"] = log_transformed.X
    rna_in.uns["log1p"] = log_transformed.uns["log1p"]
    return input_data


@pytest.fixture
def lognormed_test_data_path(tmp_path, lognormed_test_data):
    temp_h5mu = tmp_path / "lognormed.h5mu"
    lognormed_test_data.write_h5mu(temp_h5mu)
    return temp_h5mu


@pytest.fixture
def lognormed_batch_test_data_path(tmp_path, lognormed_test_data):
    temp_h5mu = tmp_path / "lognormed_batch.h5mu"
    rna_mod = lognormed_test_data.mod["rna"]
    rna_mod.obs["batch"] = "A"
    column_index = rna_mod.obs.columns.get_indexer(["batch"])
    rna_mod.obs.iloc[slice(rna_mod.n_obs // 2, None), column_index] = "B"
    lognormed_test_data.write_h5mu(temp_h5mu)
    return temp_h5mu


@pytest.fixture()
def filter_data_path(tmp_path, input_data):
    temp_h5mu = tmp_path / "filtered.h5mu"
    rna_in = input_data.mod["rna"]
    sc.pp.filter_genes(rna_in, min_counts=20)
    input_data.write_h5mu(temp_h5mu)
    return temp_h5mu


@pytest.fixture()
def common_vars_data(lognormed_test_data):
    rna_in = lognormed_test_data.mod["rna"]
    rna_in.var["common_vars"] = False
    column_index = rna_in.var.columns.get_indexer(["common_vars"])
    rna_in.var.iloc[:10000, column_index] = True
    rna_in.var["common_vars"] = rna_in.var["common_vars"].astype("boolean")
    return lognormed_test_data


@pytest.fixture()
def common_vars_data_path(tmp_path, common_vars_data):
    temp_h5mu = tmp_path / "lognormed_var_input.h5mu"
    common_vars_data.write_h5mu(temp_h5mu)
    return temp_h5mu


def test_filter_with_hvg(run_component, lognormed_test_data_path):
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            lognormed_test_data_path,
            "--output",
            "output.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
        ]
    )
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    # Put the output data back into its original shape
    # so that we can compare it to the input
    data.mod["rna"].var = data.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    data.var = data.var.drop(columns=["rna:filter_with_hvg"], errors="raise")
    del data["rna"].varm["hvg"]
    assert_annotation_objects_equal(lognormed_test_data_path, data)


def test_filter_with_hvg_var_input(run_component, common_vars_data_path):
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            common_vars_data_path,
            "--output",
            "output.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
        ]
    )

    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            common_vars_data_path,
            "--output",
            "common_vars.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
            "--var_input",
            "common_vars",
        ]
    )

    mdata = mu.read_h5mu("output.h5mu")
    common_vars = mu.read_h5mu("common_vars.h5mu")

    # Assert detected HVG are different
    hvg = mdata.mod["rna"][:, mdata.mod["rna"].var["filter_with_hvg"]].var_names
    common_vars_hvg = common_vars.mod["rna"][
        :, common_vars.mod["rna"].var["filter_with_hvg"]
    ].var_names

    assert len(hvg) != len(common_vars_hvg), (
        "Number of HVG should be different when var_input is defined"
    )

    # Assert original data is unchanged
    # Put the output data back into its original shape
    # so that we can compare it to the input
    mdata.mod["rna"].var = mdata.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    mdata.var = mdata.var.drop(columns=["rna:filter_with_hvg"], errors="raise")
    del mdata["rna"].varm["hvg"]

    common_vars.mod["rna"].var = common_vars.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    common_vars.var = common_vars.var.drop(
        columns=["rna:filter_with_hvg"], errors="raise"
    )
    del common_vars["rna"].varm["hvg"]

    assert_annotation_objects_equal(mdata, common_vars)


def test_filter_with_hvg_batch_with_batch(
    run_component, lognormed_batch_test_data_path
):
    """
    Make sure that selecting a layer works together with obs_batch_key.
    https://github.com/scverse/scanpy/issues/2396
    """
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            lognormed_batch_test_data_path,
            "--output",
            "output.h5mu",
            "--obs_batch_key",
            "batch",
            "--layer",
            "log_transformed",
        ]
    )
    assert os.path.exists("output.h5mu")
    output_data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in output_data.mod["rna"].var.columns

    # Check the contents of the output to check if the correct layer was selected
    input_mudata = mu.read_h5mu(lognormed_batch_test_data_path)
    input_data = input_mudata.mod["rna"].copy()
    input_data.X = input_data.layers["log_transformed"].copy()
    del input_data.layers["log_transformed"]
    input_data.uns["log1p"]["base"] = None
    expected_output = sc.pp.highly_variable_genes(
        input_data, batch_key="batch", inplace=False, subset=False
    )
    expected_output = expected_output.reindex(index=input_mudata.mod["rna"].var.index)
    pd.testing.assert_series_equal(
        expected_output["highly_variable"],
        output_data.mod["rna"].var["filter_with_hvg"],
        check_names=False,
    )
    # Put the output data back into its original shape
    # so that we can compare it to the input
    output_data.mod["rna"].var = output_data.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    output_data.var = output_data.var.drop(
        columns=["rna:filter_with_hvg"], errors="raise"
    )
    del output_data["rna"].varm["hvg"]
    assert_annotation_objects_equal(lognormed_batch_test_data_path, output_data)


def test_filter_with_hvg_seurat_v3_requires_n_top_features(run_component, input_path):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--flavor",
                "seurat_v3",  # Uses raw data.
                "--output",
                "output.h5mu",
            ]
        )
    assert re.search(
        "When flavor is set to 'seurat_v3', you are required to set 'n_top_features'.",
        err.value.stdout.decode("utf-8"),
    )


def test_filter_with_hvg_seurat_v3(run_component, input_path):
    run_component(
        [
            "--input",
            input_path,
            "--flavor",
            "seurat_v3",  # Uses raw data.
            "--output",
            "output.h5mu",
            "--n_top_features",
            "50",
        ]
    )
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "hvg" in data.mod["rna"].varm
    # Put the output data back into its original shape
    # so that we can compare it to the input
    data.mod["rna"].var = data.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    data.var = data.var.drop(columns=["rna:filter_with_hvg"], errors="raise")
    del data["rna"].varm["hvg"]
    assert_annotation_objects_equal(input_path, data)


def test_filter_with_hvg_cell_ranger(run_component, filter_data_path):
    run_component(
        [
            "--input",
            filter_data_path,
            "--flavor",
            "cell_ranger",  # Must use filtered data.
            "--output",
            "output.h5mu",
        ]
    )
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    # Put the output data back into its original shape
    # so that we can compare it to the input
    data.mod["rna"].var = data.mod["rna"].var.drop(
        columns=["filter_with_hvg"], errors="raise"
    )
    data.var = data.var.drop(columns=["rna:filter_with_hvg"], errors="raise")
    del data["rna"].varm["hvg"]
    assert_annotation_objects_equal(filter_data_path, data)


def test_filter_with_hvg_cell_ranger_unfiltered_data_change_error_message(
    run_component, input_path
):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--flavor",
                "cell_ranger",  # Must use filtered data, but in this test we use unfiltered data
                "--output",
                "output.h5mu",
            ]
        )
    assert re.search(
        r"Scanpy failed to calculate hvg. The error "
        r"returned by scanpy \(see above\) could be the "
        r"result from trying to use this component on unfiltered data.",
        err.value.stdout.decode("utf-8"),
    )


def test_filter_with_hvg_exclude_features(run_component, lognormed_test_data_path):
    """
    Test that excluding features from HVG calculation works correctly.
    Features in the exclusion list should be excluded from HVG selection
    but remain in the output data with highly_variable=False.
    """
    # Read input to get some feature names to exclude
    input_data = mu.read_h5mu(lognormed_test_data_path)
    rna = input_data.mod["rna"]
    # Get first 5 feature names to exclude
    features_to_exclude = rna.var_names[:5].tolist()

    # First, run without exclusion
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            lognormed_test_data_path,
            "--output",
            "output_no_exclusion.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
        ]
    )

    # Then run with exclusion of specific features
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            lognormed_test_data_path,
            "--output",
            "output_with_exclusion.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
            "--features_to_exclude",
            ";".join(features_to_exclude),
        ]
    )

    assert os.path.exists("output_with_exclusion.h5mu")
    data_no_exclusion = mu.read_h5mu("output_no_exclusion.h5mu")
    data_with_exclusion = mu.read_h5mu("output_with_exclusion.h5mu")

    rna_no_exclusion = data_no_exclusion.mod["rna"]
    rna_with_exclusion = data_with_exclusion.mod["rna"]

    # Check that filter_with_hvg column exists
    assert "filter_with_hvg" in rna_with_exclusion.var.columns

    # Check that all features are still present (exclusion doesn't remove features)
    assert rna_with_exclusion.n_vars == rna_no_exclusion.n_vars

    # Check that excluded features are not marked as HVG
    excluded_hvg = rna_with_exclusion.var.loc[features_to_exclude, "filter_with_hvg"]
    assert not excluded_hvg.any(), (
        "Excluded features should not be marked as highly variable"
    )


def test_filter_with_hvg_stores_uns(run_component, lognormed_test_data_path):
    """
    Test that the uns attribute is correctly populated with HVG information.
    """
    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            lognormed_test_data_path,
            "--output",
            "output.h5mu",
            "--layer",
            "log_transformed",
            "--output_compression",
            "gzip",
            "--uns_name",
            "hvg",
        ]
    )
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    assert "hvg" in data.mod["rna"].varm
    assert data.mod["rna"].uns["hvg"] == {"flavor": "seurat"}


def test_filter_with_hvg_var_input_with_batch(
    run_component, common_vars_data, tmp_path
):
    """
    Make sure the details stored in .varm can be written when var_input is
    combined with obs_batch_key, which adds 'highly_variable_intersection'.
    """
    input_data = common_vars_data
    rna_in = input_data.mod["rna"]
    rna_in.obs["batch"] = "A"
    column_index = rna_in.obs.columns.get_indexer(["batch"])
    rna_in.obs.iloc[slice(rna_in.n_obs // 2, None), column_index] = "B"
    input_path = tmp_path / "lognormed_var_input_batch.h5mu"
    input_data.write_h5mu(input_path)

    run_component(
        [
            "--flavor",
            "seurat",
            "--input",
            input_path,
            "--output",
            "output.h5mu",
            "--layer",
            "log_transformed",
            "--var_input",
            "common_vars",
            "--obs_batch_key",
            "batch",
        ]
    )
    assert os.path.exists("output.h5mu")
    data = mu.read_h5mu("output.h5mu")
    hvg = data.mod["rna"].varm["hvg"]
    common_vars = data.mod["rna"].var["common_vars"].to_numpy(dtype=bool)
    for column in ("highly_variable", "highly_variable_intersection"):
        assert hvg[column].dtype == bool
        assert not hvg.loc[~common_vars, column].any()


@pytest.mark.parametrize("flavor", ["seurat", "cell_ranger"])
@pytest.mark.parametrize("subset_by", ["var_input", "features_to_exclude"])
def test_filter_with_hvg_subset_with_batch_matches_scanpy(
    run_component, common_vars_data, tmp_path, flavor, subset_by
):
    """
    With obs_batch_key, scanpy returns the dispersion based flavors sorted by
    feature name instead of in input order. Make sure the results are still
    assigned to the correct features when only a subset of the features is used.
    """
    input_data = common_vars_data
    rna_in = input_data.mod["rna"]
    rna_in.obs["batch"] = "A"
    column_index = rna_in.obs.columns.get_indexer(["batch"])
    rna_in.obs.iloc[slice(rna_in.n_obs // 2, None), column_index] = "B"
    input_path = tmp_path / "lognormed_var_input_batch.h5mu"
    input_data.write_h5mu(input_path)

    if subset_by == "var_input":
        keep = rna_in.var["common_vars"].to_numpy(dtype=bool)
        subset_args = ["--var_input", "common_vars"]
    else:
        # Leave out few enough features to keep --features_to_exclude within
        # the command line length limit.
        keep = (np.arange(rna_in.n_vars) % 50) != 0
        subset_args = ["--features_to_exclude", ";".join(rna_in.var_names[~keep])]

    run_component(
        [
            "--flavor",
            flavor,
            "--input",
            input_path,
            "--output",
            "output.h5mu",
            "--layer",
            "log_transformed",
            "--obs_batch_key",
            "batch",
            *subset_args,
        ]
    )
    assert os.path.exists("output.h5mu")
    output_rna = mu.read_h5mu("output.h5mu").mod["rna"]

    # Calculate the expected result by running scanpy on the subset directly
    expected_input = ad.AnnData(
        X=rna_in.layers["log_transformed"].copy(),
        obs=rna_in.obs[["batch"]],
        var=pd.DataFrame(index=rna_in.var_names),
        uns={"log1p": {"base": None}},
    )[:, keep].copy()
    expected = sc.pp.highly_variable_genes(
        expected_input, flavor=flavor, batch_key="batch", inplace=False, subset=False
    )
    kept_features = rna_in.var_names[keep]
    expected = expected.reindex(index=kept_features)
    assert expected["highly_variable"].any()

    hvg = output_rna.varm["hvg"]
    hvg.index = output_rna.var_names
    for column in ("highly_variable", "dispersions_norm", "highly_variable_nbatches"):
        pd.testing.assert_series_equal(
            hvg.loc[kept_features, column],
            expected[column],
            check_names=False,
            check_dtype=False,
        )
    pd.testing.assert_series_equal(
        output_rna.var.loc[kept_features, "filter_with_hvg"],
        expected["highly_variable"],
        check_names=False,
    )
    assert not output_rna.var.loc[~keep, "filter_with_hvg"].any()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
