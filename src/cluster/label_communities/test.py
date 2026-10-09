import sys
import subprocess
import pytest
import numpy as np
import pandas as pd

## VIASH START
meta = {
    "executable": "target/executable/cluster/label_communities/label_communities",
    "resources_dir": "src/cluster/label_communities/",
    "config": "src/cluster/label_communities/config.vsh.yaml",
}
## VIASH END

ID = "participant_id"


def _make_tables(tmp_path, n_groups=60, seed=42, id_column=ID):
    """Proportion and dynamics tables with two planted communities.

    Labels `a1`/`a2`/`a3` move together across groups; `b1`/`b2`/`b3` move
    together in the opposite direction. Any working community detection must
    recover exactly that split, so a broken similarity matrix cannot pass.
    """
    rng = np.random.default_rng(seed)
    ids = [f"donor_{i:03d}" for i in range(n_groups)]
    pt = np.linspace(0.0, 1.0, n_groups)

    group_a = np.column_stack([0.2 + 0.5 * pt + rng.normal(0, 0.01, n_groups)] * 3)
    group_b = np.column_stack([0.7 - 0.5 * pt + rng.normal(0, 0.01, n_groups)] * 3)
    # De-correlate the members within a community slightly
    group_a += rng.normal(0, 0.005, group_a.shape)
    group_b += rng.normal(0, 0.005, group_b.shape)

    labels = ["a1", "a2", "a3", "b1", "b2", "b3"]
    raw = np.clip(np.hstack([group_a, group_b]), 1e-6, None)
    props = raw / raw.sum(axis=1, keepdims=True)

    prop_df = pd.DataFrame(props, columns=labels)
    prop_df.insert(0, id_column, ids)
    prop_path = tmp_path / "proportions.csv"
    prop_df.to_csv(prop_path, index=False)

    # Fitted curves: community a rises, community b falls
    grid = np.linspace(0.0, 1.0, 30)
    rows = []
    for label in labels:
        curve = 0.2 + 0.5 * grid if label.startswith("a") else 0.7 - 0.5 * grid
        for x, y in zip(grid, curve):
            rows.append({"label": label, "pseudotime": x, "proportion_fitted": y})
    dyn_path = tmp_path / "dynamics.csv"
    pd.DataFrame(rows).to_csv(dyn_path, index=False)

    return prop_path, dyn_path, labels


def _communities(path):
    df = pd.read_csv(path)
    return dict(zip(df["label"], df["community_id"].astype(str)))


def test_basic(run_component, tmp_path):
    """One row per label, with the planted split recovered."""
    prop_path, dyn_path, labels = _make_tables(tmp_path)
    output = tmp_path / "communities.csv"

    run_component(
        [
            "--input",
            str(prop_path),
            "--dynamics",
            str(dyn_path),
            "--n_communities",
            "2",
            "--output",
            str(output),
        ]
    )

    assert output.is_file()
    result = pd.read_csv(output)
    assert list(result.columns) == ["label", "community_id"]
    assert sorted(result["label"]) == sorted(labels)

    comm = _communities(output)
    assert comm["a1"] == comm["a2"] == comm["a3"], f"Community a split: {comm}"
    assert comm["b1"] == comm["b2"] == comm["b3"], f"Community b split: {comm}"
    assert comm["a1"] != comm["b1"], f"Communities a and b merged: {comm}"


@pytest.mark.parametrize("method", ["hierarchical", "spectral"])
def test_methods(run_component, tmp_path, method):
    """Both clustering methods recover the planted split."""
    prop_path, dyn_path, _ = _make_tables(tmp_path)
    output = tmp_path / f"communities_{method}.csv"

    run_component(
        [
            "--input",
            str(prop_path),
            "--dynamics",
            str(dyn_path),
            "--n_communities",
            "2",
            "--method",
            method,
            "--output",
            str(output),
        ]
    )

    comm = _communities(output)
    assert comm["a1"] == comm["a3"] != comm["b1"], f"{method} failed: {comm}"


@pytest.mark.parametrize("correlation", ["pearson", "spearman"])
def test_correlation_method(run_component, tmp_path, correlation):
    """Both correlation methods are accepted and recover the split."""
    prop_path, dyn_path, _ = _make_tables(tmp_path)
    output = tmp_path / f"communities_{correlation}.csv"

    run_component(
        [
            "--input",
            str(prop_path),
            "--dynamics",
            str(dyn_path),
            "--n_communities",
            "2",
            "--correlation_method",
            correlation,
            "--output",
            str(output),
        ]
    )

    comm = _communities(output)
    assert comm["a2"] != comm["b2"], f"{correlation} failed: {comm}"


def test_alpha_one_without_dynamics(run_component, tmp_path):
    """--alpha 1.0 is co-occurrence only, so --dynamics may be omitted."""
    prop_path, _, _ = _make_tables(tmp_path)
    output = tmp_path / "communities_alpha1.csv"

    run_component(
        [
            "--input",
            str(prop_path),
            "--n_communities",
            "2",
            "--alpha",
            "1.0",
            "--output",
            str(output),
        ]
    )

    comm = _communities(output)
    assert comm["a1"] == comm["a2"] != comm["b1"], f"Co-occurrence only failed: {comm}"


def test_output_similarity(run_component, tmp_path):
    """--output_similarity holds one row per label pair."""
    prop_path, dyn_path, labels = _make_tables(tmp_path)
    output = tmp_path / "communities.csv"
    sim_path = tmp_path / "similarity.csv"

    run_component(
        [
            "--input",
            str(prop_path),
            "--dynamics",
            str(dyn_path),
            "--n_communities",
            "2",
            "--output",
            str(output),
            "--output_similarity",
            str(sim_path),
        ]
    )

    sim = pd.read_csv(sim_path)
    n = len(labels)
    assert list(sim.columns) == [
        "label_1",
        "label_2",
        "co_occurrence",
        "dynamics",
        "combined",
    ]
    assert len(sim) == n * (n - 1) // 2
    assert sim[["co_occurrence", "dynamics", "combined"]].abs().max().max() <= 1.0

    within = sim[
        sim["label_1"].str.startswith("a") & sim["label_2"].str.startswith("a")
    ]["co_occurrence"].mean()
    between = sim[
        sim["label_1"].str.startswith("a") & sim["label_2"].str.startswith("b")
    ]["co_occurrence"].mean()
    assert within > between, (
        f"Within-community co-occurrence {within:.3f} not above between {between:.3f}"
    )


def test_missing_dynamics(run_component, tmp_path):
    """--dynamics is required whenever --alpha is below 1.0."""
    prop_path, _, _ = _make_tables(tmp_path)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(prop_path),
                "--n_communities",
                "2",
                "--output",
                str(tmp_path / "out.csv"),
            ]
        )
    assert "--dynamics is required" in err.value.stdout.decode("utf-8")


def test_too_many_communities(run_component, tmp_path):
    """--n_communities above the number of labels is rejected."""
    prop_path, dyn_path, labels = _make_tables(tmp_path)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(prop_path),
                "--dynamics",
                str(dyn_path),
                "--n_communities",
                str(len(labels) + 1),
                "--output",
                str(tmp_path / "out.csv"),
            ]
        )
    assert "exceeds the number of labels" in err.value.stdout.decode("utf-8")


def test_dynamics_missing_column(run_component, tmp_path):
    """A --dynamics table without the expected columns is rejected."""
    prop_path, dyn_path, _ = _make_tables(tmp_path)
    bad = pd.read_csv(dyn_path).rename(columns={"proportion_fitted": "value"})
    bad_path = tmp_path / "dynamics_bad.csv"
    bad.to_csv(bad_path, index=False)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(prop_path),
                "--dynamics",
                str(bad_path),
                "--n_communities",
                "2",
                "--output",
                str(tmp_path / "out.csv"),
            ]
        )
    assert "missing column" in err.value.stdout.decode("utf-8")


def test_zero_variance_label(run_component, tmp_path):
    """A label with constant proportions gets zero similarity, not a crash."""
    prop_path, dyn_path, _ = _make_tables(tmp_path)
    df = pd.read_csv(prop_path)
    df["constant"] = 0.05
    padded = tmp_path / "proportions_constant.csv"
    df.to_csv(padded, index=False)

    output = tmp_path / "communities_constant.csv"
    run_component(
        [
            "--input",
            str(padded),
            "--dynamics",
            str(dyn_path),
            "--n_communities",
            "3",
            "--output",
            str(output),
        ]
    )

    result = pd.read_csv(output)
    assert "constant" in set(result["label"])
    assert not result["community_id"].isna().any()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
