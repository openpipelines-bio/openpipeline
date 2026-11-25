import sys
import subprocess
import pytest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/trajectory/fit_proportion_dynamics/fit_proportion_dynamics",
    "resources_dir": "src/trajectory/fit_proportion_dynamics/",
    "config": "src/trajectory/fit_proportion_dynamics/config.vsh.yaml",
}
## VIASH END

LABELS = ["rising", "falling", "peaking"]


def _make_mudata(
    tmp_path, n_groups=40, seed=42, pseudotime_column="palantir_pseudotime"
):
    """Group-level MuData with known proportion dynamics.

    Three labels with different shapes along pseudotime, so a broken fit is
    visible rather than merely noisy:
      - `rising`  monotonically increases
      - `falling` monotonically decreases
      - `peaking` peaks in the middle
    Proportions of the three sum to 1 per group, as real compositions do.
    """
    rng = np.random.default_rng(seed)
    ids = [f"donor_{i:03d}" for i in range(n_groups)]
    pt = np.linspace(0.0, 1.0, n_groups)

    rising = 0.2 + 0.6 * pt
    falling = 0.8 - 0.6 * pt
    peaking = 0.4 - 1.2 * (pt - 0.5) ** 2
    raw = np.column_stack([rising, falling, peaking])
    raw = raw + rng.normal(0, 0.01, raw.shape)
    raw = np.clip(raw, 1e-6, None)
    props = raw / raw.sum(axis=1, keepdims=True)

    obs = pd.DataFrame(
        {pseudotime_column: pt, "palantir_entropy": rng.random(n_groups)}, index=ids
    )
    adata = ad.AnnData(X=props, obs=obs, var=pd.DataFrame(index=LABELS))
    # Shuffled, so the component cannot rely on the groups being ordered by pseudotime
    adata = adata[rng.permutation(n_groups)].copy()

    path = tmp_path / "proportions.h5mu"
    mu.MuData({"proportions": adata}).write_h5mu(str(path))
    return path


def _read(path):
    return mu.read_h5mu(str(path)).mod["proportions"]


def test_basic(run_component, tmp_path):
    """One fitted curve per label, on the requested number of bins."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "dynamics.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--output",
            str(output),
            "--n_pseudotime_bins",
            "50",
        ]
    )

    assert output.is_file()
    adata = _read(output)
    assert list(adata.var_names) == LABELS
    curves = adata.varm["proportion_dynamics"]
    assert curves.shape == (3, 50)
    assert np.isfinite(curves).all()
    grid = adata.uns["proportion_dynamics_pseudotime"]
    assert len(grid) == 50
    assert np.all(np.diff(grid) > 0), "pseudotime bins are not increasing"
    # the proportions themselves are left alone
    np.testing.assert_allclose(adata.X, _read(input_path).X)


def test_curve_shapes(run_component, tmp_path):
    """The fitted curves and statistics reproduce the shapes that were put in."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "dynamics.h5mu"

    run_component(
        ["--input", str(input_path), "--output", str(output), "--lam", "0.001"]
    )

    adata = _read(output)
    curves = pd.DataFrame(adata.varm["proportion_dynamics"], index=adata.var_names)
    assert curves.loc["rising"].iloc[-1] > curves.loc["rising"].iloc[0], (
        "'rising' does not rise"
    )
    assert curves.loc["falling"].iloc[-1] < curves.loc["falling"].iloc[0], (
        "'falling' does not fall"
    )

    stats = adata.var
    assert [c for c in stats.columns if c.startswith("dynamics_")] == [
        "dynamics_peak_pseudotime",
        "dynamics_r_squared",
        "dynamics_p_value",
        "dynamics_n_groups",
    ]
    assert stats.loc["peaking", "dynamics_peak_pseudotime"] == pytest.approx(
        0.5, abs=0.2
    ), f"'peaking' peaks at {stats.loc['peaking', 'dynamics_peak_pseudotime']}"
    assert stats.loc["rising", "dynamics_peak_pseudotime"] > 0.8
    assert stats.loc["falling", "dynamics_peak_pseudotime"] < 0.2
    assert (stats["dynamics_n_groups"] == 40).all()


def test_missing_pseudotime(run_component, tmp_path):
    """Groups without a pseudotime value are left out of the fit."""
    input_path = _make_mudata(tmp_path)
    mdata = mu.read_h5mu(str(input_path))
    mdata.mod["proportions"].obs.iloc[0, 0] = np.nan
    trimmed = tmp_path / "proportions_trimmed.h5mu"
    mdata.write_h5mu(str(trimmed))

    output = tmp_path / "dynamics.h5mu"
    run_component(["--input", str(trimmed), "--output", str(output)])

    n_groups = _read(output).var["dynamics_n_groups"]
    assert (n_groups == 39).all(), (
        f"Expected 39 groups after dropping one, got {n_groups.unique()}"
    )


def test_custom_keys(run_component, tmp_path):
    """--obs_pseudotime, --varm_output, --uns_output and --var_prefix are respected."""
    input_path = _make_mudata(tmp_path, pseudotime_column="via_time")
    output = tmp_path / "dynamics.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obs_pseudotime",
            "via_time",
            "--varm_output",
            "my_curves",
            "--uns_output",
            "my_grid",
            "--var_prefix",
            "my_",
            "--output",
            str(output),
        ]
    )

    adata = _read(output)
    assert adata.varm["my_curves"].shape == (3, 100)
    assert len(adata.uns["my_grid"]) == 100
    assert "my_p_value" in adata.var.columns


@pytest.mark.parametrize(
    "mutation,message",
    [
        ("bad_pseudotime_column", "--obs_pseudotime 'nope' not found"),
        ("non_numeric", "must be numeric"),
        ("too_few_groups", "--min_groups"),
    ],
)
def test_input_errors(run_component, tmp_path, mutation, message):
    """Each malformed input is rejected with a specific message."""
    input_path = _make_mudata(tmp_path)
    args = ["--input", str(input_path), "--output", str(tmp_path / "out.h5mu")]

    if mutation == "bad_pseudotime_column":
        args += ["--obs_pseudotime", "nope"]
    elif mutation == "non_numeric":
        args += ["--obs_pseudotime", "label"]
        mdata = mu.read_h5mu(str(input_path))
        mdata.mod["proportions"].obs["label"] = "a"
        mdata.write_h5mu(str(input_path))
    elif mutation == "too_few_groups":
        mdata = mu.read_h5mu(str(input_path))
        mdata.mod["proportions"].obs.iloc[3:, 0] = np.nan
        mdata.write_h5mu(str(input_path))

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert message in err.value.stdout.decode("utf-8")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
