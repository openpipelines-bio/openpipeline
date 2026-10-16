from __future__ import annotations
import sys
import subprocess
import pytest
import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu

## VIASH START
meta = {
    "executable": "./target/executable/metadata/calculate_obs_label_proportions/calculate_obs_label_proportions",
    "resources_dir": "./src/metadata/calculate_obs_label_proportions/",
    "cpus": 2,
    "config": "./src/metadata/calculate_obs_label_proportions/config.vsh.yaml",
}
## VIASH END


def _make_mudata(tmp_path, n_donors=3, n_subpops=4, cells_per_group=10):
    """Create a synthetic MuData object with participant_id and subpopulation columns."""
    participants = [f"donor_{i}" for i in range(n_donors)]
    subpopulations = [f"subpop_{j}" for j in range(n_subpops)]

    obs_rows = []
    for pid in participants:
        for subpop in subpopulations:
            for _ in range(cells_per_group):
                obs_rows.append({"participant_id": pid, "subpopulation": subpop})

    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell_{k}" for k in range(len(obs))]

    # No X: the component only reads .obs
    adata = AnnData(obs=obs)
    mdata = MuData({"rna": adata})

    h5mu_path = tmp_path / "input.h5mu"
    mdata.write_h5mu(str(h5mu_path))
    return mdata, h5mu_path


@pytest.fixture
def synthetic_mudata(tmp_path):
    return _make_mudata(tmp_path)


def _read_proportions(path, modality="proportions"):
    """Proportion matrix of a group-level output, as a groups x labels DataFrame."""
    adata = read_h5mu(str(path)).mod[modality]
    return pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)


def test_basic(run_component, tmp_path):
    """One observation per group and one variable per label, in a new modality."""
    mdata, h5mu_path = _make_mudata(
        tmp_path, n_donors=3, n_subpops=4, cells_per_group=10
    )
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obs_group",
            "participant_id",
            "--obs_label",
            "subpopulation",
        ]
    )

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    assert list(result.mod) == ["proportions"], (
        f"expected a single 'proportions' modality, got {list(result.mod)}"
    )
    adata = result.mod["proportions"]

    # one observation per group, one variable per label
    assert adata.shape == (3, 4)
    assert sorted(adata.obs_names) == ["donor_0", "donor_1", "donor_2"]
    assert sorted(adata.var_names) == [f"subpop_{j}" for j in range(4)]
    assert list(adata.obs["participant_id"]) == list(adata.obs_names)

    # Each row must sum to 1
    np.testing.assert_allclose(adata.X.sum(axis=1), np.ones(3), atol=1e-10)

    # obsm copy for components that read .obsm, and the raw cell counts
    np.testing.assert_array_equal(adata.obsm["X_proportions"], adata.X)
    np.testing.assert_array_equal(adata.layers["counts"], np.full((3, 4), 10))


def test_uniform_proportions(run_component, tmp_path):
    """With equal cell counts per group, all proportions should be 1/n_subpops."""
    n_subpops = 4
    mdata, h5mu_path = _make_mudata(
        tmp_path, n_donors=2, n_subpops=n_subpops, cells_per_group=5
    )
    output_path = tmp_path / "output_uniform.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obs_group",
            "participant_id",
            "--obs_label",
            "subpopulation",
        ]
    )

    prop_df = _read_proportions(output_path)
    expected = 1.0 / n_subpops
    np.testing.assert_allclose(prop_df.values, expected, atol=1e-10)


def test_custom_column_names(run_component, tmp_path):
    """Component works with non-default obs column names."""
    obs = pd.DataFrame(
        {
            "sample": ["s1"] * 10 + ["s2"] * 10,
            "cell_class": (["A"] * 5 + ["B"] * 5) * 2,
        },
        index=[f"c{i}" for i in range(20)],
    )
    adata = AnnData(obs=obs)
    mdata = MuData({"rna": adata})
    h5mu_path = tmp_path / "custom_cols.h5mu"
    mdata.write_h5mu(str(h5mu_path))
    output_path = tmp_path / "output_custom.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obs_group",
            "sample",
            "--obs_label",
            "cell_class",
            "--output_modality",
            "my_props",
            "--obsm_output",
            "my_obsm",
        ]
    )

    result = read_h5mu(str(output_path))
    assert list(result.mod) == ["my_props"]
    adata_out = result.mod["my_props"]
    assert "my_obsm" in adata_out.obsm
    assert list(adata_out.obs["sample"]) == ["s1", "s2"]
    np.testing.assert_allclose(adata_out.X.sum(axis=1), np.ones(2), atol=1e-10)


def test_normalize_within(run_component, tmp_path):
    """--obs_normalize_within makes each label's denominator its own class.

    Two classes of very different size. Globally the small class's labels are a few
    percent of the donor; within their own class they are half of it. The two
    normalisations therefore cannot be confused for one another.
    """
    composition = {
        "donor_A": {
            ("big", "b1"): 90,
            ("big", "b2"): 10,
            ("small", "s1"): 5,
            ("small", "s2"): 5,
        },
        "donor_B": {
            ("big", "b1"): 25,
            ("big", "b2"): 75,
            ("small", "s1"): 8,
            ("small", "s2"): 2,
        },
    }
    obs_rows = []
    for donor, cells in composition.items():
        for (cls, label), n in cells.items():
            obs_rows += [
                {"participant_id": donor, "cell_class": cls, "subpopulation": label}
            ] * n
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell_{k}" for k in range(len(obs))]
    h5mu_path = tmp_path / "hierarchy.h5mu"
    MuData({"rna": AnnData(obs=obs)}).write_h5mu(str(h5mu_path))
    output_path = tmp_path / "output_within.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obs_group",
            "participant_id",
            "--obs_label",
            "subpopulation",
            "--obs_normalize_within",
            "cell_class",
        ]
    )

    prop = _read_proportions(output_path)
    expected = {
        "donor_A": {"b1": 0.9, "b2": 0.1, "s1": 0.5, "s2": 0.5},
        "donor_B": {"b1": 0.25, "b2": 0.75, "s1": 0.8, "s2": 0.2},
    }
    for donor, values in expected.items():
        for label, value in values.items():
            assert prop.loc[donor, label] == pytest.approx(value, abs=1e-12), (
                f"{donor}/{label}: expected {value}, got {prop.loc[donor, label]}"
            )

    # Two classes, so a row sums to 2 and not to 1
    np.testing.assert_allclose(prop.sum(axis=1).to_numpy(), np.full(2, 2.0), atol=1e-12)

    # The label -> class hierarchy is kept in .var
    var = read_h5mu(str(output_path)).mod["proportions"].var
    assert var["cell_class"].to_dict() == {
        "b1": "big",
        "b2": "big",
        "s1": "small",
        "s2": "small",
    }


def test_normalize_within_differs_from_global(run_component, tmp_path):
    """Without --obs_normalize_within the same input gives the global proportions."""
    obs_rows = [
        {"participant_id": "d", "cell_class": "big", "subpopulation": "b1"}
    ] * 90 + [
        {"participant_id": "d", "cell_class": "small", "subpopulation": "s1"}
    ] * 10
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell_{k}" for k in range(len(obs))]
    h5mu_path = tmp_path / "global.h5mu"
    MuData({"rna": AnnData(obs=obs)}).write_h5mu(str(h5mu_path))
    output_path = tmp_path / "output_global.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obs_group",
            "participant_id",
            "--obs_label",
            "subpopulation",
        ]
    )

    prop = _read_proportions(output_path)
    assert prop.loc["d", "b1"] == pytest.approx(0.9)
    assert prop.loc["d", "s1"] == pytest.approx(0.1)


def test_normalize_within_requires_hierarchy(run_component, tmp_path):
    """A label that spans two classes is rejected, not silently double-counted."""
    obs_rows = [
        {"participant_id": "d", "cell_class": "big", "subpopulation": "shared"}
    ] * 10 + [
        {"participant_id": "d", "cell_class": "small", "subpopulation": "shared"}
    ] * 10
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell_{k}" for k in range(len(obs))]
    h5mu_path = tmp_path / "ambiguous.h5mu"
    MuData({"rna": AnnData(obs=obs)}).write_h5mu(str(h5mu_path))

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(h5mu_path),
                "--output",
                str(tmp_path / "out.h5mu"),
                "--obs_group",
                "participant_id",
                "--obs_label",
                "subpopulation",
                "--obs_normalize_within",
                "cell_class",
            ]
        )
    assert "expects a strict hierarchy" in err.value.stdout.decode("utf-8")


def test_missing_column_raises(run_component, tmp_path):
    """An informative error is raised when the obs column is missing."""
    mdata, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_err.h5mu"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(h5mu_path),
                "--output",
                str(output_path),
                "--obs_group",
                "participant_id",
                "--obs_label",
                "nonexistent_column",
            ]
        )
    assert "Column 'nonexistent_column' not found in .obs" in err.value.stdout.decode(
        "utf-8"
    )


def test_output_compression(run_component, tmp_path):
    """Component writes compressed output without error."""
    mdata, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_compressed.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obs_group",
            "participant_id",
            "--obs_label",
            "subpopulation",
            "--output_compression",
            "gzip",
        ]
    )

    assert output_path.is_file()
    assert _read_proportions(output_path).shape == (3, 4)


def test_unequal_group_sizes(run_component, tmp_path):
    """Proportions must be right when donors have different cell counts per subpopulation.

    Equal-sized groups make every proportion 1/n_subpops, which a wrong denominator
    would also produce; the composition below is asymmetric in both directions and a
    subpopulation is absent from one donor.
    """
    composition = {
        "donor_A": {"sub_0": 10, "sub_1": 30, "sub_2": 0},
        "donor_B": {"sub_0": 5, "sub_1": 5, "sub_2": 10},
        "donor_C": {"sub_0": 1, "sub_1": 0, "sub_2": 3},
    }
    expected = {
        "donor_A": {"sub_0": 0.25, "sub_1": 0.75, "sub_2": 0.0},
        "donor_B": {"sub_0": 0.25, "sub_1": 0.25, "sub_2": 0.5},
        "donor_C": {"sub_0": 0.25, "sub_1": 0.0, "sub_2": 0.75},
    }

    obs_rows = []
    for donor, subpops in composition.items():
        for subpop, n in subpops.items():
            obs_rows += [{"participant_id": donor, "subpopulation": subpop}] * n
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell_{k}" for k in range(len(obs))]

    h5mu_path = tmp_path / "unequal.h5mu"
    MuData({"rna": AnnData(obs=obs)}).write_h5mu(str(h5mu_path))
    output_path = tmp_path / "output_unequal.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obs_group",
            "participant_id",
            "--obs_label",
            "subpopulation",
        ]
    )

    prop_df = _read_proportions(output_path)
    expected_df = pd.DataFrame(expected).T.loc[prop_df.index, prop_df.columns]

    np.testing.assert_allclose(
        prop_df.to_numpy(dtype=float),
        expected_df.to_numpy(dtype=float),
        atol=1e-12,
        err_msg=f"Proportions differ from hand-computed values:\n{prop_df}",
    )


def test_overall_proportions_without_group(run_component, tmp_path):
    """Without --obs_group all cells form one group named 'all'."""
    obs_rows = (
        [{"participant_id": "donor_A", "subpopulation": "sub_0"}] * 10
        + [{"participant_id": "donor_B", "subpopulation": "sub_0"}] * 20
        + [{"participant_id": "donor_B", "subpopulation": "sub_1"}] * 10
    )
    obs = pd.DataFrame(obs_rows, index=[f"cell_{k}" for k in range(40)])
    h5mu_path = tmp_path / "overall.h5mu"
    MuData({"rna": AnnData(obs=obs)}).write_h5mu(str(h5mu_path))
    output_path = tmp_path / "output_overall.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obs_label",
            "subpopulation",
        ]
    )

    adata_out = read_h5mu(str(output_path)).mod["proportions"]
    prop_df = _read_proportions(output_path)
    assert list(prop_df.index) == ["all"], f"expected a single 'all' row, got {prop_df}"
    # 30/40 cells are sub_0, 10/40 are sub_1 - not the per-donor proportions
    np.testing.assert_allclose(prop_df.loc["all", "sub_0"], 0.75, atol=1e-12)
    np.testing.assert_allclose(prop_df.loc["all", "sub_1"], 0.25, atol=1e-12)
    # the helper column must not leak into the written object
    assert list(adata_out.obs.columns) == ["group"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
