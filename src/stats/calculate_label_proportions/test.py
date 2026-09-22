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
    "executable": "./target/executable/stats/calculate_label_proportions/calculate_label_proportions",
    "resources_dir": "./src/stats/calculate_label_proportions/",
    "cpus": 2,
    "config": "./src/stats/calculate_label_proportions/config.vsh.yaml",
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


def test_basic(run_component, tmp_path):
    """Proportions are computed and stored in uns, and in obsm when asked for."""
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
            "--obsm_output",
            "proportions",
        ]
    )

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    adata = result.mod["rna"]

    # uns["proportions"] is a DataFrame (groups x labels), not a nested dict
    assert "proportions" in adata.uns
    prop_df = adata.uns["proportions"]
    assert isinstance(prop_df, pd.DataFrame), (
        f"uns['proportions'] should survive the round-trip as a DataFrame, got {type(prop_df)}"
    )
    assert prop_df.shape == (3, 4)

    # Each row must sum to 1
    row_sums = prop_df.sum(axis=1)
    np.testing.assert_allclose(row_sums.values, np.ones(3), atol=1e-10)

    # obsm["proportions"] must exist with correct shape
    assert "proportions" in adata.obsm
    obsm = adata.obsm["proportions"]
    assert obsm.shape == (len(adata), 4)


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

    result = read_h5mu(str(output_path))
    prop_dict = result.mod["rna"].uns["proportions"]
    prop_df = pd.DataFrame(prop_dict)
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
            "--uns_output",
            "my_props",
            "--obsm_output",
            "my_props",
        ]
    )

    result = read_h5mu(str(output_path))
    adata_out = result.mod["rna"]
    assert "my_props" in adata_out.uns
    assert "my_props" in adata_out.obsm
    prop_df = pd.DataFrame(adata_out.uns["my_props"]).T
    np.testing.assert_allclose(prop_df.sum(axis=1).values, np.ones(2), atol=1e-10)


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
    result = read_h5mu(str(output_path))
    assert "proportions" in result.mod["rna"].uns


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
            "--obsm_output",
            "proportions",
        ]
    )

    result = read_h5mu(str(output_path))
    adata_out = result.mod["rna"]
    prop_df = pd.DataFrame(adata_out.uns["proportions"])
    expected_df = pd.DataFrame(expected).T.loc[prop_df.index, prop_df.columns]

    np.testing.assert_allclose(
        prop_df.to_numpy(dtype=float),
        expected_df.to_numpy(dtype=float),
        atol=1e-12,
        err_msg=f"Proportions differ from hand-computed values:\n{prop_df}",
    )

    # every cell carries its own donor's row
    obsm = pd.DataFrame(
        np.asarray(adata_out.obsm["proportions"]),
        index=adata_out.obs_names,
        columns=prop_df.columns,
    )
    for donor in composition:
        rows = obsm[adata_out.obs["participant_id"].to_numpy() == donor]
        np.testing.assert_allclose(
            rows.to_numpy(dtype=float),
            np.tile(
                expected_df.loc[donor].to_numpy(dtype=float), (len(rows), 1)
            ),
            atol=1e-12,
            err_msg=f"obsm rows for {donor} do not match its proportion vector",
        )


def test_obsm_not_written_by_default(run_component, tmp_path):
    """The redundant per-cell copy is only written when --obsm_output is given."""
    mdata, h5mu_path = _make_mudata(tmp_path, n_donors=2, n_subpops=3, cells_per_group=5)
    output_path = tmp_path / "output_no_obsm.h5mu"

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

    adata_out = read_h5mu(str(output_path)).mod["rna"]
    assert "proportions" in adata_out.uns
    assert "proportions" not in adata_out.obsm, (
        f".obsm should be empty without --obsm_output, found {list(adata_out.obsm)}"
    )


def test_output_csv(run_component, tmp_path):
    """--output_csv writes the same matrix as a table, keyed by the group column."""
    mdata, h5mu_path = _make_mudata(
        tmp_path, n_donors=3, n_subpops=4, cells_per_group=10
    )
    output_path = tmp_path / "output_csv.h5mu"
    csv_path = tmp_path / "proportions.csv"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--output_csv",
            str(csv_path),
            "--obs_group",
            "participant_id",
            "--obs_label",
            "subpopulation",
        ]
    )

    assert csv_path.is_file()
    csv_df = pd.read_csv(str(csv_path), index_col=0)
    assert csv_df.index.name == "participant_id"
    assert csv_df.shape == (3, 4)
    np.testing.assert_allclose(csv_df.sum(axis=1).values, np.ones(3), atol=1e-10)

    # identical to what went into .uns
    uns_df = pd.DataFrame(read_h5mu(str(output_path)).mod["rna"].uns["proportions"])
    np.testing.assert_allclose(
        csv_df.loc[uns_df.index, uns_df.columns].to_numpy(dtype=float),
        uns_df.to_numpy(dtype=float),
        atol=1e-12,
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

    adata_out = read_h5mu(str(output_path)).mod["rna"]
    prop_df = pd.DataFrame(adata_out.uns["proportions"])
    assert list(prop_df.index) == ["all"], f"expected a single 'all' row, got {prop_df}"
    # 30/40 cells are sub_0, 10/40 are sub_1 - not the per-donor proportions
    np.testing.assert_allclose(prop_df.loc["all", "sub_0"], 0.75, atol=1e-12)
    np.testing.assert_allclose(prop_df.loc["all", "sub_1"], 0.25, atol=1e-12)
    # the helper column must not leak into the written object
    assert "_obs_group" not in adata_out.obs.columns


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
