import sys
import subprocess
import pytest
import numpy as np
import pandas as pd

## VIASH START
meta = {
    "executable": "target/executable/trajectory/fit_proportion_dynamics/fit_proportion_dynamics",
    "resources_dir": "src/trajectory/fit_proportion_dynamics/",
    "config": "src/trajectory/fit_proportion_dynamics/config.vsh.yaml",
}
## VIASH END

ID = "participant_id"


def _make_tables(tmp_path, n_groups=40, seed=42, id_column=ID):
    """Two group-level tables with known proportion dynamics.

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

    prop_df = pd.DataFrame(props, columns=["rising", "falling", "peaking"])
    prop_df.insert(0, id_column, ids)
    prop_path = tmp_path / "proportions.csv"
    prop_df.to_csv(prop_path, index=False)

    # Shuffled, so the component cannot rely on row order to align the tables
    pt_df = pd.DataFrame(
        {
            id_column: ids,
            "palantir_pseudotime": pt,
            "palantir_entropy": rng.random(n_groups),
        }
    ).sample(frac=1.0, random_state=seed)
    pt_path = tmp_path / "pseudotime.csv"
    pt_df.to_csv(pt_path, index=False)

    return prop_path, pt_path, ids


def test_basic(run_component, tmp_path):
    """Long-format curves for every label, on the requested number of bins."""
    prop_path, pt_path, _ = _make_tables(tmp_path)
    output = tmp_path / "dynamics.csv"

    run_component(
        [
            "--input",
            str(prop_path),
            "--pseudotime",
            str(pt_path),
            "--output",
            str(output),
            "--n_pseudotime_bins",
            "50",
        ]
    )

    assert output.is_file()
    result = pd.read_csv(output)
    assert list(result.columns) == ["label", "pseudotime", "proportion_fitted"]
    assert sorted(result["label"].unique()) == ["falling", "peaking", "rising"]
    assert len(result) == 3 * 50
    assert np.isfinite(result["proportion_fitted"]).all()


def test_curve_shapes(run_component, tmp_path):
    """The fitted curves reproduce the shapes that were put in."""
    prop_path, pt_path, _ = _make_tables(tmp_path)
    output = tmp_path / "dynamics.csv"
    stats_path = tmp_path / "stats.csv"

    run_component(
        [
            "--input",
            str(prop_path),
            "--pseudotime",
            str(pt_path),
            "--output",
            str(output),
            "--output_stats",
            str(stats_path),
            "--lam",
            "0.001",
        ]
    )

    curves = pd.read_csv(output)
    by_label = {
        label: group.sort_values("pseudotime")["proportion_fitted"].to_numpy()
        for label, group in curves.groupby("label")
    }
    assert by_label["rising"][-1] > by_label["rising"][0], "'rising' does not rise"
    assert by_label["falling"][-1] < by_label["falling"][0], "'falling' does not fall"

    stats = pd.read_csv(stats_path).set_index("label")
    assert list(stats.columns) == [
        "peak_pseudotime",
        "r_squared",
        "p_value",
        "n_groups",
    ]
    assert stats.loc["peaking", "peak_pseudotime"] == pytest.approx(0.5, abs=0.2), (
        f"'peaking' peaks at {stats.loc['peaking', 'peak_pseudotime']}, expected ~0.5"
    )
    assert stats.loc["rising", "peak_pseudotime"] > 0.8
    assert stats.loc["falling", "peak_pseudotime"] < 0.2
    assert (stats["n_groups"] == 40).all()


def test_partial_overlap(run_component, tmp_path):
    """Groups present in only one table are dropped, the rest still fit."""
    prop_path, pt_path, ids = _make_tables(tmp_path)
    trimmed = pd.read_csv(pt_path)
    trimmed = trimmed[trimmed[ID] != ids[0]]
    trimmed_path = tmp_path / "pseudotime_trimmed.csv"
    trimmed.to_csv(trimmed_path, index=False)

    output = tmp_path / "dynamics.csv"
    stats_path = tmp_path / "stats.csv"
    run_component(
        [
            "--input",
            str(prop_path),
            "--pseudotime",
            str(trimmed_path),
            "--output",
            str(output),
            "--output_stats",
            str(stats_path),
        ]
    )

    stats = pd.read_csv(stats_path)
    assert (stats["n_groups"] == 39).all(), (
        f"Expected 39 groups after dropping one, got {stats['n_groups'].unique()}"
    )


def test_custom_pseudotime_column(run_component, tmp_path):
    """--pseudotime_column selects a differently named column."""
    prop_path, pt_path, _ = _make_tables(tmp_path)
    renamed = pd.read_csv(pt_path).rename(columns={"palantir_pseudotime": "via_time"})
    renamed_path = tmp_path / "pseudotime_renamed.csv"
    renamed.to_csv(renamed_path, index=False)

    output = tmp_path / "dynamics.csv"
    run_component(
        [
            "--input",
            str(prop_path),
            "--pseudotime",
            str(renamed_path),
            "--pseudotime_column",
            "via_time",
            "--output",
            str(output),
        ]
    )
    assert pd.read_csv(output)["label"].nunique() == 3


def test_id_column(run_component, tmp_path):
    """--id_column picks a non-first identifier column shared by both tables."""
    prop_path, pt_path, _ = _make_tables(tmp_path, id_column="donor")
    prop = pd.read_csv(prop_path)
    prop = prop[[c for c in prop.columns if c != "donor"] + ["donor"]]
    reordered = tmp_path / "proportions_reordered.csv"
    prop.to_csv(reordered, index=False)

    output = tmp_path / "dynamics.csv"
    run_component(
        [
            "--input",
            str(reordered),
            "--pseudotime",
            str(pt_path),
            "--id_column",
            "donor",
            "--output",
            str(output),
        ]
    )
    assert pd.read_csv(output)["label"].nunique() == 3


@pytest.mark.parametrize(
    "mutation,message",
    [
        ("duplicate_id", "duplicated identifiers"),
        ("non_numeric", "non-numeric value column"),
        ("bad_id_column", "not found in --pseudotime"),
        ("bad_pseudotime_column", "--pseudotime_column"),
        ("too_few_groups", "--min_groups"),
    ],
)
def test_input_errors(run_component, tmp_path, mutation, message):
    """Each malformed input is rejected with a specific message."""
    prop_path, pt_path, _ = _make_tables(tmp_path)
    args = [
        "--input",
        str(prop_path),
        "--pseudotime",
        str(pt_path),
        "--output",
        str(tmp_path / "out.csv"),
    ]

    if mutation == "duplicate_id":
        df = pd.read_csv(prop_path)
        df.loc[1, ID] = df.loc[0, ID]
        path = tmp_path / "dup.csv"
        df.to_csv(path, index=False)
        args[1] = str(path)
    elif mutation == "non_numeric":
        df = pd.read_csv(prop_path)
        df["batch"] = "a"
        path = tmp_path / "non_numeric.csv"
        df.to_csv(path, index=False)
        args[1] = str(path)
    elif mutation == "bad_id_column":
        df = pd.read_csv(pt_path).rename(columns={ID: "other_id"})
        path = tmp_path / "bad_id.csv"
        df.to_csv(path, index=False)
        args[3] = str(path)
    elif mutation == "bad_pseudotime_column":
        args += ["--pseudotime_column", "nope"]
    elif mutation == "too_few_groups":
        df = pd.read_csv(pt_path).head(3)
        path = tmp_path / "short.csv"
        df.to_csv(path, index=False)
        args[3] = str(path)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert message in err.value.stdout.decode("utf-8")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
