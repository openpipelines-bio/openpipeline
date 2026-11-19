import sys
import subprocess
import pytest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/trajectory/palantir/palantir",
    "resources_dir": "src/trajectory/palantir/",
    "config": "src/trajectory/palantir/config.vsh.yaml",
}
## VIASH END


def _make_mudata(tmp_path, n_cells=300, n_dims=10, seed=42):
    """Synthetic MuData with a linear trajectory in X_pca_integrated.

    Cells are arranged along a 1-D pseudotime axis embedded in n_dims dimensions
    with small noise, giving Palantir a clear gradient to recover.
    """
    rng = np.random.default_rng(seed)
    pt = np.linspace(0, 1, n_cells)

    # Linear trajectory + noise across all dims
    embedding = np.column_stack(
        [pt]
        + [
            pt * rng.uniform(0.5, 1.0) + rng.normal(0, 0.05, n_cells)
            for _ in range(n_dims - 1)
        ]
    ).astype("float32")

    # No X: Palantir only reads the embedding in .obsm
    cluster = np.repeat(["A", "B", "C", "D", "E"], n_cells // 5)[:n_cells]
    obs = pd.DataFrame(
        {"cluster": cluster},
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )

    adata = ad.AnnData(obs=obs)
    adata.obsm["X_pca_integrated"] = embedding

    mdata = mu.MuData({"rna": adata})
    path = tmp_path / "palantir_input.h5mu"
    mdata.write_h5mu(str(path))
    return path


def test_palantir_start_id_barcode(run_component, tmp_path):
    """Run Palantir using an explicit --start_id barcode."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obsm_input",
            "X_pca_integrated",
            "--start_id",
            "cell_0000",
            "--num_waypoints",
            "50",
            "--n_components",
            "5",
            "--knn",
            "10",
            "--output",
            str(output),
        ]
    )

    assert output.is_file(), "Output h5mu not created"

    adata = mu.read_h5mu(str(output)).mod["rna"]

    assert "palantir_pseudotime" in adata.obs.columns
    assert "palantir_entropy" in adata.obs.columns
    assert "palantir_fate_probabilities" in adata.obsm
    assert "palantir_waypoints" in adata.uns

    pt = adata.obs["palantir_pseudotime"]
    assert pt.notna().all(), "Pseudotime contains NaN values"
    assert (pt >= 0).all() and (pt <= 1).all(), "Pseudotime outside [0, 1]"
    assert len(adata.uns["palantir_waypoints"]) > 0, "Waypoints list is empty"


def test_palantir_start_cluster(run_component, tmp_path):
    """Start cell is resolved automatically from a cluster label."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "output_cluster.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obsm_input",
            "X_pca_integrated",
            "--start_cluster",
            "A",
            "--start_cluster_column",
            "cluster",
            "--num_waypoints",
            "50",
            "--n_components",
            "5",
            "--knn",
            "10",
            "--output",
            str(output),
        ]
    )

    assert output.is_file()
    adata = mu.read_h5mu(str(output)).mod["rna"]
    assert "palantir_pseudotime" in adata.obs.columns
    assert adata.obs["palantir_pseudotime"].notna().all()


def test_palantir_custom_output_keys(run_component, tmp_path):
    """Custom obs/obsm/uns output key names are respected."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "output_custom.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obsm_input",
            "X_pca_integrated",
            "--start_id",
            "cell_0000",
            "--num_waypoints",
            "50",
            "--n_components",
            "5",
            "--knn",
            "10",
            "--pseudotime_column",
            "my_pseudotime",
            "--entropy_column",
            "my_entropy",
            "--obsm_fate_probabilities",
            "my_fate_probs",
            "--uns_waypoints",
            "my_waypoints",
            "--output",
            str(output),
        ]
    )

    adata = mu.read_h5mu(str(output)).mod["rna"]
    assert "my_pseudotime" in adata.obs.columns
    assert "my_entropy" in adata.obs.columns
    assert "my_fate_probs" in adata.obsm
    assert "my_waypoints" in adata.uns


def test_palantir_input_preserved(run_component, tmp_path):
    """Cell barcodes and gene names are unchanged after Palantir."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "output_preserved.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obsm_input",
            "X_pca_integrated",
            "--start_id",
            "cell_0000",
            "--num_waypoints",
            "50",
            "--n_components",
            "5",
            "--knn",
            "10",
            "--output",
            str(output),
        ]
    )

    orig = mu.read_h5mu(str(input_path)).mod["rna"]
    out = mu.read_h5mu(str(output)).mod["rna"]

    np.testing.assert_array_equal(orig.obs_names, out.obs_names)
    np.testing.assert_array_equal(orig.var_names, out.var_names)


# ---------------------------------------------------------------------------
# Group-level input: one observation per donor, as written by
# metadata/calculate_obs_label_proportions + dimred/phate
# ---------------------------------------------------------------------------


def _make_group_mudata(tmp_path, n_groups=120, n_dims=4, seed=42):
    """Group-level MuData with a linear trajectory in X_phate and a stage label."""
    rng = np.random.default_rng(seed)
    pt = np.linspace(0, 1, n_groups)
    embedding = np.column_stack(
        [pt]
        + [
            pt * rng.uniform(0.5, 1.0) + rng.normal(0, 0.05, n_groups)
            for _ in range(n_dims - 1)
        ]
    )
    ids = [f"donor_{i:03d}" for i in range(n_groups)]

    # Early third of the trajectory is the "control" stage
    stage = np.where(pt < 0.33, "control", np.where(pt < 0.66, "mid", "late"))
    obs = pd.DataFrame({"participant_id": ids, "stage": stage}, index=ids)
    adata = ad.AnnData(obs=obs)
    adata.obsm["X_phate"] = embedding

    path = tmp_path / "groups.h5mu"
    mu.MuData({"proportions": adata}).write_h5mu(str(path))
    return path, ids


def _group_args(input_path, output, *extra):
    return [
        "--input",
        str(input_path),
        "--modality",
        "proportions",
        "--obsm_input",
        "X_phate",
        "--output",
        str(output),
        "--num_waypoints",
        "40",
        "--waypoint_knn",
        "10",
        "--n_components",
        "3",
        "--knn",
        "15",
        *extra,
    ]


def test_group_level_start_id(run_component, tmp_path):
    """Explicit --start_id on donors; pseudotime follows the gradient."""
    input_path, ids = _make_group_mudata(tmp_path)
    output = tmp_path / "groups_pt.h5mu"

    run_component(_group_args(input_path, output, "--start_id", ids[0]))

    adata = mu.read_h5mu(str(output)).mod["proportions"]
    assert list(adata.obs_names) == ids, "Identifiers not preserved"
    pt = adata.obs["palantir_pseudotime"].to_numpy()
    assert np.isfinite(pt).all()
    assert pt.min() >= 0.0 and pt.max() <= 1.0
    # The donors are ordered along the trajectory, so pseudotime must increase
    rank_corr = pd.Series(pt).corr(pd.Series(np.arange(len(pt))), method="spearman")
    assert rank_corr > 0.9, (
        f"Pseudotime does not follow the trajectory: rho={rank_corr:.3f}"
    )


def test_group_level_start_cluster(run_component, tmp_path):
    """--start_cluster resolves against an .obs column of the group-level object."""
    input_path, ids = _make_group_mudata(tmp_path)
    output = tmp_path / "groups_cluster.h5mu"

    run_component(
        _group_args(
            input_path,
            output,
            "--start_cluster_column",
            "stage",
            "--start_cluster",
            "control",
        )
    )

    pt = mu.read_h5mu(str(output)).mod["proportions"].obs["palantir_pseudotime"]
    n_control = int(len(ids) * 0.33)
    assert pt.iloc[:n_control].mean() < pt.iloc[-n_control:].mean(), (
        "Root was not taken from the 'control' stage"
    )


def test_group_level_terminal_states(run_component, tmp_path):
    """Explicit --terminal_states give one fate-probability column each."""
    input_path, ids = _make_group_mudata(tmp_path)
    output = tmp_path / "groups_fates.h5mu"

    run_component(
        _group_args(
            input_path,
            output,
            "--start_id",
            ids[0],
            "--terminal_states",
            ids[-1],
            "--terminal_states",
            ids[-2],
        )
    )

    fates = (
        mu.read_h5mu(str(output)).mod["proportions"].obsm["palantir_fate_probabilities"]
    )
    assert sorted(fates.columns) == sorted([ids[-1], ids[-2]]), (
        f"Unexpected fate columns: {list(fates.columns)}"
    )
    np.testing.assert_allclose(fates.to_numpy().sum(axis=1), 1.0, atol=1e-6)


def test_waypoint_knn_too_large(run_component, tmp_path):
    """--waypoint_knn above the waypoint count fails with a message naming it."""
    input_path, ids = _make_group_mudata(tmp_path, n_groups=40)
    args = _group_args(input_path, tmp_path / "out.h5mu", "--start_id", ids[0])
    args[args.index("--num_waypoints") + 1] = "12"
    args[args.index("--waypoint_knn") + 1] = "20"
    args[args.index("--knn") + 1] = "10"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    stdout = err.value.stdout.decode("utf-8")
    assert "--waypoint_knn (20) must be smaller than the number of waypoints" in stdout


def test_waypoints_capped_at_n_obs(run_component, tmp_path):
    """--num_waypoints above the observation count is capped, and the cap is checked."""
    input_path, ids = _make_group_mudata(tmp_path, n_groups=25)
    args = _group_args(input_path, tmp_path / "out.h5mu", "--start_id", ids[0])
    args[args.index("--num_waypoints") + 1] = "500"
    args[args.index("--waypoint_knn") + 1] = "30"
    args[args.index("--knn") + 1] = "10"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    stdout = err.value.stdout.decode("utf-8")
    assert "25 = min(--num_waypoints 500, 25 observations)" in stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
