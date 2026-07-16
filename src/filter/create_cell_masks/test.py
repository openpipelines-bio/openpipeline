import re
import sys
from subprocess import CalledProcessError

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest

## VIASH START
meta = {
    "executable": "./target/executable/filter/create_cell_masks/create_cell_masks",
    "resources_dir": "src/utils",
    "config": "./src/filter/create_cell_masks/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def input_h5mu():
    rng = np.random.default_rng(seed=1)
    n_obs = 100

    rna = ad.AnnData(
        X=rng.integers(0, 100, size=(n_obs, 10)),
        obs=pd.DataFrame(
            {
                "total_counts": rng.integers(0, 100, size=n_obs),
                "control_probe_counts": rng.integers(0, 20, size=n_obs),
                "preexisting": rng.random(n_obs),
            },
            index=[f"cell_{i}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(
            {"gene_meta": [f"g_{i}" for i in range(10)]},
            index=[f"gene_{i}" for i in range(10)],
        ),
    )
    rna.obsm["existing_obsm"] = rng.random((n_obs, 3))
    rna.uns["existing_uns"] = {"params": 42}

    # A second modality to confirm cross-modality untouched.
    prot = ad.AnnData(
        X=rng.integers(0, 100, size=(n_obs, 5)),
        obs=pd.DataFrame(index=rna.obs_names),
        var=pd.DataFrame(index=[f"prot_{i}" for i in range(5)]),
    )

    return mu.MuData({"rna": rna, "prot": prot})


@pytest.fixture
def input_path(input_h5mu, random_h5mu_path):
    path = random_h5mu_path()
    input_h5mu.write(path)
    return path


def _assert_untouched(
    adata_before, adata_after, *, new_obs_cols, new_obsm_keys, new_uns_keys
):
    """Confirm that only the expected new keys appear; nothing else was changed."""
    np.testing.assert_array_equal(adata_before.X, adata_after.X)
    pd.testing.assert_index_equal(adata_before.obs_names, adata_after.obs_names)
    pd.testing.assert_index_equal(adata_before.var_names, adata_after.var_names)
    pd.testing.assert_frame_equal(adata_before.var, adata_after.var)
    # Pre-existing obs columns survive identically.
    for col in adata_before.obs.columns:
        pd.testing.assert_series_equal(
            adata_before.obs[col], adata_after.obs[col], check_names=False
        )
    # Only the new obs columns are added; nothing else.
    added_obs = set(adata_after.obs.columns) - set(adata_before.obs.columns)
    assert added_obs == set(new_obs_cols), f"unexpected obs cols: {added_obs}"
    # Pre-existing obsm survives identically.
    for key in adata_before.obsm:
        np.testing.assert_array_equal(adata_before.obsm[key], adata_after.obsm[key])
    added_obsm = set(adata_after.obsm.keys()) - set(adata_before.obsm.keys())
    assert added_obsm == set(new_obsm_keys), f"unexpected obsm keys: {added_obsm}"
    # Pre-existing uns survives.
    for key in adata_before.uns:
        assert key in adata_after.uns
    added_uns = set(adata_after.uns.keys()) - set(adata_before.uns.keys())
    assert added_uns == set(new_uns_keys), f"unexpected uns keys: {added_uns}"


def test_with_prefix_and_groups(
    run_component, input_h5mu, input_path, random_h5mu_path
):
    output = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output,
            "--modality",
            "rna",
            "--filters",
            "total_counts:gt:10:rna;control_probe_counts:le:10:control;total_counts:le:50:rna",
            "--prefix",
            "cell",
        ]
    )
    assert output.is_file()

    mdata = mu.read_h5mu(str(output))
    rna = mdata["rna"]

    # Individual masks DataFrame in .obsm (prefixed).
    masks_df = rna.obsm["cell_masks"]
    assert isinstance(masks_df, pd.DataFrame)
    assert len(masks_df.columns) == 6  # 3 individual + 2 group + overall
    for name in (
        "rna_total_counts_gt_10",
        "control_control_probe_counts_le_10",
        "rna_total_counts_le_50",
        "rna",
        "control",
        "overall",
    ):
        assert name in masks_df.columns

    # Group masks in .obs.
    for col in ("cell_mask", "cell_mask_rna", "cell_mask_control"):
        assert col in rna.obs.columns
        assert rna.obs[col].dtype == bool

    # Overall mask is the AND of all individual masks.
    expected_overall = (
        masks_df["rna_total_counts_gt_10"]
        & masks_df["control_control_probe_counts_le_10"]
        & masks_df["rna_total_counts_le_50"]
    )
    np.testing.assert_array_equal(
        rna.obs["cell_mask"].to_numpy(), expected_overall.to_numpy()
    )

    # Filter definitions in .uns.
    filters_df = rna.uns["cell_filters"]
    assert isinstance(filters_df, pd.DataFrame)
    assert len(filters_df) == 3
    assert set(filters_df["group"].tolist()) == {"rna", "control"}

    # Nothing else on the RNA modality changed.
    _assert_untouched(
        input_h5mu["rna"],
        rna,
        new_obs_cols={"cell_mask", "cell_mask_rna", "cell_mask_control"},
        new_obsm_keys={"cell_masks"},
        new_uns_keys={"cell_filters"},
    )

    # The other modality is untouched.
    prot_before = input_h5mu["prot"]
    prot_after = mdata["prot"]
    np.testing.assert_array_equal(prot_before.X, prot_after.X)
    pd.testing.assert_frame_equal(prot_before.obs, prot_after.obs)
    assert set(prot_after.obsm.keys()) == set(prot_before.obsm.keys())


def test_filter_without_group(run_component, input_h5mu, input_path, random_h5mu_path):
    """A filter without a group still produces an individual mask and contributes
    to the 'overall' group mask, but produces no per-group mask column."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output,
            "--modality",
            "rna",
            "--filters",
            "total_counts:gt:20",
            "--prefix",
            "cell",
        ]
    )

    rna = mu.read_h5mu(str(output))["rna"]
    masks_df = rna.obsm["cell_masks"]

    # Individual mask exists, no per-group, just 'overall'.
    assert "total_counts_gt_20" in masks_df.columns
    assert "overall" in masks_df.columns
    # No spurious group columns (only the two we expect).
    assert set(masks_df.columns) == {"total_counts_gt_20", "overall"}

    # No group-specific obs column.
    obs_mask_cols = [c for c in rna.obs.columns if c.startswith("cell_mask")]
    assert obs_mask_cols == ["cell_mask"]

    filters_df = rna.uns["cell_filters"]
    assert filters_df["group"].iloc[0] == ""


def test_without_prefix(run_component, input_h5mu, input_path, random_h5mu_path):
    """Omitting --prefix yields un-prefixed slot names."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output,
            "--modality",
            "rna",
            "--filters",
            "total_counts:gt:10:rna;control_probe_counts:le:10:control",
        ]
    )

    rna = mu.read_h5mu(str(output))["rna"]
    assert "masks" in rna.obsm
    assert "filters" in rna.uns
    for col in ("mask", "mask_rna", "mask_control"):
        assert col in rna.obs.columns

    _assert_untouched(
        input_h5mu["rna"],
        rna,
        new_obs_cols={"mask", "mask_rna", "mask_control"},
        new_obsm_keys={"masks"},
        new_uns_keys={"filters"},
    )


def test_missing_column_errors_fast(run_component, input_path, random_h5mu_path):
    """Column existence is validated up-front for all filters."""
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--modality",
                "rna",
                "--filters",
                "total_counts:gt:10;not_a_column:gt:5",
            ]
        )
    assert re.search(r"not_a_column", err.value.stdout.decode("utf-8"))


def test_bad_operator_errors(run_component, input_path, random_h5mu_path):
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                random_h5mu_path(),
                "--modality",
                "rna",
                "--filters",
                "total_counts:foo:10",
            ]
        )
    assert re.search(r"Unknown operator 'foo'", err.value.stdout.decode("utf-8"))


def test_compression(run_component, input_path, random_h5mu_path):
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input_path,
            "--output",
            output,
            "--modality",
            "rna",
            "--filters",
            "total_counts:gt:10:rna",
            "--output_compression",
            "gzip",
        ]
    )
    assert output.is_file()
    rna = mu.read_h5mu(str(output))["rna"]
    assert "masks" in rna.obsm


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
