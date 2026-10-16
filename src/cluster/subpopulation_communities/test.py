import sys
import subprocess
import pytest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/cluster/subpopulation_communities/subpopulation_communities",
    "resources_dir": "src/cluster/subpopulation_communities/",
    "config": "src/cluster/subpopulation_communities/config.vsh.yaml",
}
## VIASH END

LABELS = ["a1", "a2", "a3", "b1", "b2", "b3"]


def _make_mudata(tmp_path, n_groups=60, seed=42, with_dynamics=True, extra=None):
    """Group-level MuData with two planted communities.

    Labels `a1`/`a2`/`a3` move together across groups; `b1`/`b2`/`b3` move
    together in the opposite direction. Any working community detection must
    recover exactly that split, so a broken similarity matrix cannot pass.
    `extra` maps additional label names to (proportions, fitted curve).
    """
    rng = np.random.default_rng(seed)
    ids = [f"donor_{i:03d}" for i in range(n_groups)]
    pt = np.linspace(0.0, 1.0, n_groups)

    group_a = np.column_stack([0.2 + 0.5 * pt + rng.normal(0, 0.01, n_groups)] * 3)
    group_b = np.column_stack([0.7 - 0.5 * pt + rng.normal(0, 0.01, n_groups)] * 3)
    # De-correlate the members within a community slightly
    group_a += rng.normal(0, 0.005, group_a.shape)
    group_b += rng.normal(0, 0.005, group_b.shape)

    raw = np.clip(np.hstack([group_a, group_b]), 1e-6, None)
    props = raw / raw.sum(axis=1, keepdims=True)

    # Fitted curves: community a rises, community b falls
    grid = np.linspace(0.0, 1.0, 30)
    curves = np.vstack(
        [
            0.2 + 0.5 * grid if label.startswith("a") else 0.7 - 0.5 * grid
            for label in LABELS
        ]
    )

    labels = list(LABELS)
    for label, (values, curve) in (extra or {}).items():
        labels.append(label)
        props = np.column_stack([props, values])
        curves = np.vstack([curves, curve])

    adata = ad.AnnData(
        X=props, obs=pd.DataFrame(index=ids), var=pd.DataFrame(index=labels)
    )
    if with_dynamics:
        adata.varm["proportion_dynamics"] = curves
        adata.uns["proportion_dynamics_pseudotime"] = grid

    path = tmp_path / "proportions.h5mu"
    mu.MuData({"proportions": adata}).write_h5mu(str(path))
    return path, labels


def _communities(path):
    var = mu.read_h5mu(str(path)).mod["proportions"].var
    return var["community_id"].astype(str).to_dict()


def test_basic(run_component, tmp_path):
    """One community per label, with the planted split recovered."""
    input_path, labels = _make_mudata(tmp_path)
    output = tmp_path / "communities.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--n_communities",
            "2",
            "--output",
            str(output),
        ]
    )

    assert output.is_file()
    comm = _communities(output)
    assert sorted(comm) == sorted(labels)
    assert comm["a1"] == comm["a2"] == comm["a3"], f"Community a split: {comm}"
    assert comm["b1"] == comm["b2"] == comm["b3"], f"Community b split: {comm}"
    assert comm["a1"] != comm["b1"], f"Communities a and b merged: {comm}"


@pytest.mark.parametrize("method", ["hierarchical", "spectral"])
def test_methods(run_component, tmp_path, method):
    """Both clustering methods recover the planted split."""
    input_path, _ = _make_mudata(tmp_path)
    output = tmp_path / f"communities_{method}.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
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
    input_path, _ = _make_mudata(tmp_path)
    output = tmp_path / f"communities_{correlation}.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
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
    """--alpha 1.0 is co-occurrence only, so no fitted curves are needed."""
    input_path, _ = _make_mudata(tmp_path, with_dynamics=False)
    output = tmp_path / "communities_alpha1.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
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


def test_similarity_matrices(run_component, tmp_path):
    """The three label x label similarity matrices are written to .varp."""
    input_path, labels = _make_mudata(tmp_path)
    output = tmp_path / "communities.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--n_communities",
            "2",
            "--output",
            str(output),
        ]
    )

    adata = mu.read_h5mu(str(output)).mod["proportions"]
    n = len(labels)
    for key in ("co_occurrence", "dynamics", "combined"):
        sim = np.asarray(adata.varp[f"community_{key}"])
        assert sim.shape == (n, n), f"{key}: shape {sim.shape}"
        np.testing.assert_allclose(sim, sim.T, atol=1e-12)
        assert np.abs(sim).max() <= 1.0

    co = pd.DataFrame(
        np.asarray(adata.varp["community_co_occurrence"]), index=labels, columns=labels
    )
    within = co.loc[["a1", "a2", "a3"], ["a1", "a2", "a3"]].to_numpy()
    within = within[~np.eye(3, dtype=bool)].mean()
    between = co.loc[["a1", "a2", "a3"], ["b1", "b2", "b3"]].to_numpy().mean()
    assert within > between, (
        f"Within-community co-occurrence {within:.3f} not above between {between:.3f}"
    )


def test_missing_dynamics(run_component, tmp_path):
    """The fitted curves are required whenever --alpha is below 1.0."""
    input_path, _ = _make_mudata(tmp_path, with_dynamics=False)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_path),
                "--n_communities",
                "2",
                "--output",
                str(tmp_path / "out.h5mu"),
            ]
        )
    assert "--varm_dynamics 'proportion_dynamics' not found" in (
        err.value.stdout.decode("utf-8")
    )


def test_too_many_communities(run_component, tmp_path):
    """--n_communities above the number of labels is rejected."""
    input_path, labels = _make_mudata(tmp_path)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_path),
                "--n_communities",
                str(len(labels) + 1),
                "--output",
                str(tmp_path / "out.h5mu"),
            ]
        )
    assert "exceeds the number of labels" in err.value.stdout.decode("utf-8")


def test_unfitted_label(run_component, tmp_path):
    """A label whose curve could not be fitted (NaN) gets zero dynamics similarity."""
    rng = np.random.default_rng(3)
    input_path, _ = _make_mudata(
        tmp_path,
        extra={"unfitted": (rng.uniform(0.01, 0.02, 60), np.full(30, np.nan))},
    )
    output = tmp_path / "communities_unfitted.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--n_communities",
            "3",
            "--output",
            str(output),
        ]
    )

    adata = mu.read_h5mu(str(output)).mod["proportions"]
    dyn = pd.DataFrame(
        np.asarray(adata.varp["community_dynamics"]),
        index=adata.var_names,
        columns=adata.var_names,
    )
    assert (dyn.loc["unfitted"].drop("unfitted") == 0).all()
    assert not adata.var["community_id"].isna().any()


def test_zero_variance_label(run_component, tmp_path):
    """A label with constant proportions gets zero similarity, not a crash."""
    input_path, _ = _make_mudata(
        tmp_path, extra={"constant": (np.full(60, 0.05), np.full(30, 0.05))}
    )
    output = tmp_path / "communities_constant.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--n_communities",
            "3",
            "--output",
            str(output),
        ]
    )

    comm = _communities(output)
    assert "constant" in comm
    assert "nan" not in comm.values()


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
