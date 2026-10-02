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
# Table mode: group-level input, no MuData involved
# ---------------------------------------------------------------------------


def _make_tables(tmp_path, n_groups=120, n_dims=4, seed=42, id_column="participant_id"):
    """Embedding CSV with a linear trajectory, plus a metadata CSV of labels."""
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

    emb_df = pd.DataFrame(embedding, columns=[f"phate_{i + 1}" for i in range(n_dims)])
    emb_df.insert(0, id_column, ids)
    emb_path = tmp_path / "embedding.csv"
    emb_df.to_csv(emb_path, index=False)

    # Early third of the trajectory is the "control" stage
    stage = np.where(pt < 0.33, "control", np.where(pt < 0.66, "mid", "late"))
    meta_df = pd.DataFrame({id_column: ids, "stage": stage, "age": pt * 40 + 50})
    meta_path = tmp_path / "metadata.csv"
    meta_df.to_csv(meta_path, index=False)

    return emb_path, meta_path, ids


def test_table_mode_start_id(run_component, tmp_path):
    """Explicit --start_id identifier on a table; pseudotime follows the gradient."""
    emb_path, _, ids = _make_tables(tmp_path)
    output = tmp_path / "pseudotime.csv"

    run_component(
        [
            "--input_table",
            str(emb_path),
            "--output_table",
            str(output),
            "--start_id",
            ids[0],
            "--num_waypoints",
            "40",
            "--waypoint_knn",
            "10",
            "--n_components",
            "3",
            "--knn",
            "15",
        ]
    )

    assert output.is_file()
    result = pd.read_csv(output)
    assert result.columns[0] == "participant_id"
    assert "palantir_pseudotime" in result.columns
    assert "palantir_entropy" in result.columns
    assert "palantir_waypoint" in result.columns
    assert list(result["participant_id"]) == ids, "Identifiers not preserved"

    pt = result["palantir_pseudotime"].to_numpy()
    assert np.isfinite(pt).all()
    assert pt.min() >= 0.0 and pt.max() <= 1.0
    # The table rows are ordered along the trajectory, so pseudotime must increase
    rank_corr = pd.Series(pt).corr(pd.Series(np.arange(len(pt))), method="spearman")
    assert rank_corr > 0.9, (
        f"Pseudotime does not follow the trajectory: rho={rank_corr:.3f}"
    )


def test_table_mode_start_cluster_from_metadata(run_component, tmp_path):
    """--start_cluster resolves against a column of --metadata."""
    emb_path, meta_path, ids = _make_tables(tmp_path)
    output = tmp_path / "pseudotime_cluster.csv"

    run_component(
        [
            "--input_table",
            str(emb_path),
            "--metadata",
            str(meta_path),
            "--output_table",
            str(output),
            "--start_cluster_column",
            "stage",
            "--start_cluster",
            "control",
            "--num_waypoints",
            "40",
            "--waypoint_knn",
            "10",
            "--n_components",
            "3",
            "--knn",
            "15",
        ]
    )

    result = pd.read_csv(output)
    pt = result.set_index("participant_id")["palantir_pseudotime"]
    n_control = int(len(ids) * 0.33)
    assert pt.iloc[:n_control].mean() < pt.iloc[-n_control:].mean(), (
        "Root was not taken from the 'control' stage"
    )


def test_table_mode_terminal_states(run_component, tmp_path):
    """Explicit --terminal_states produce one fate_<id> column each."""
    emb_path, _, ids = _make_tables(tmp_path)
    output = tmp_path / "pseudotime_fates.csv"

    run_component(
        [
            "--input_table",
            str(emb_path),
            "--output_table",
            str(output),
            "--start_id",
            ids[0],
            "--terminal_states",
            ids[-1],
            "--terminal_states",
            ids[-2],
            "--num_waypoints",
            "40",
            "--waypoint_knn",
            "10",
            "--n_components",
            "3",
            "--knn",
            "15",
        ]
    )

    result = pd.read_csv(output)
    fate_cols = [c for c in result.columns if c.startswith("fate_")]
    assert sorted(fate_cols) == sorted([f"fate_{ids[-1]}", f"fate_{ids[-2]}"]), (
        f"Unexpected fate columns: {fate_cols}"
    )
    fates = result[fate_cols].to_numpy()
    np.testing.assert_allclose(fates.sum(axis=1), 1.0, atol=1e-6)


def test_table_mode_metadata_missing_id(run_component, tmp_path):
    """--metadata that does not cover every identifier is an error."""
    emb_path, meta_path, _ = _make_tables(tmp_path, n_groups=40)
    truncated = pd.read_csv(meta_path).head(10)
    truncated_path = tmp_path / "metadata_short.csv"
    truncated.to_csv(truncated_path, index=False)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_table",
                str(emb_path),
                "--metadata",
                str(truncated_path),
                "--output_table",
                str(tmp_path / "out.csv"),
                "--start_cluster_column",
                "stage",
                "--start_cluster",
                "control",
                "--num_waypoints",
                "15",
                "--waypoint_knn",
                "5",
                "--n_components",
                "3",
                "--knn",
                "10",
            ]
        )
    assert "is missing" in err.value.stdout.decode("utf-8")


def test_metadata_rejected_with_h5mu(run_component, tmp_path):
    """--metadata only applies to table mode."""
    input_path = _make_mudata(tmp_path, n_cells=60, n_dims=5)
    _, meta_path, _ = _make_tables(tmp_path, n_groups=10)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_path),
                "--metadata",
                str(meta_path),
                "--output",
                str(tmp_path / "out.h5mu"),
                "--start_id",
                "cell_0000",
            ]
        )
    assert "--metadata only applies to --input_table" in err.value.stdout.decode(
        "utf-8"
    )


@pytest.mark.parametrize(
    "args,message",
    [
        (
            ["--input", "H5MU", "--input_table", "TABLE", "--output_table", "OUT"],
            "Exactly one of",
        ),
        ([], "Exactly one of"),
        (["--input_table", "TABLE"], "--output_table is required"),
        (["--input", "H5MU"], "--output is required"),
    ],
)
def test_input_mode_errors(run_component, tmp_path, args, message):
    """Exactly one input mode, with the matching output argument."""
    input_path = _make_mudata(tmp_path, n_cells=60, n_dims=5)
    emb_path, _, _ = _make_tables(tmp_path, n_groups=20)
    substitutions = {
        "TABLE": str(emb_path),
        "H5MU": str(input_path),
        "OUT": str(tmp_path / "out.csv"),
    }
    resolved = [substitutions.get(a, a) for a in args]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(resolved)
    assert message in err.value.stdout.decode("utf-8")


def test_waypoint_knn_too_large(run_component, tmp_path):
    """--waypoint_knn above the waypoint count fails with a message naming it."""
    emb_path, _, ids = _make_tables(tmp_path, n_groups=40)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_table",
                str(emb_path),
                "--output_table",
                str(tmp_path / "out.csv"),
                "--start_id",
                ids[0],
                "--num_waypoints",
                "12",
                "--waypoint_knn",
                "20",
                "--n_components",
                "3",
                "--knn",
                "10",
            ]
        )
    stdout = err.value.stdout.decode("utf-8")
    assert "--waypoint_knn (20) must be smaller than the number of waypoints" in stdout


def test_waypoints_capped_at_n_obs(run_component, tmp_path):
    """--num_waypoints above the row count is capped, and the cap is what is checked."""
    emb_path, _, ids = _make_tables(tmp_path, n_groups=25)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input_table",
                str(emb_path),
                "--output_table",
                str(tmp_path / "out.csv"),
                "--start_id",
                ids[0],
                "--num_waypoints",
                "500",
                "--waypoint_knn",
                "30",
                "--n_components",
                "3",
                "--knn",
                "10",
            ]
        )
    stdout = err.value.stdout.decode("utf-8")
    assert "25 = min(--num_waypoints 500, 25 observations)" in stdout


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
