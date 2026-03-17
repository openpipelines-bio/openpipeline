from __future__ import annotations
import sys
import pytest
import numpy as np
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu

## VIASH START
meta = {
    "executable": "./target/executable/metadata/calculate_proportions/calculate_proportions",
    "resources_dir": "./src/metadata/calculate_proportions/",
    "cpus": 2,
    "config": "./src/metadata/calculate_proportions/config.vsh.yaml",
}
## VIASH END


def _make_mudata(tmp_path, n_donors=3, n_subpops=4, cells_per_group=10):
    """Create a synthetic MuData object with participant_id and subpopulation columns."""
    rng = np.random.default_rng(42)
    participants = [f"donor_{i}" for i in range(n_donors)]
    subpopulations = [f"subpop_{j}" for j in range(n_subpops)]

    obs_rows = []
    for pid in participants:
        for subpop in subpopulations:
            for _ in range(cells_per_group):
                obs_rows.append({"participant_id": pid, "subpopulation": subpop})

    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell_{k}" for k in range(len(obs))]

    n_cells = len(obs)
    n_genes = 20
    X = rng.integers(0, 100, size=(n_cells, n_genes)).astype(float)
    var = pd.DataFrame(index=[f"gene_{g}" for g in range(n_genes)])

    adata = AnnData(X=X, obs=obs, var=var)
    mdata = MuData({"rna": adata})

    h5mu_path = tmp_path / "input.h5mu"
    mdata.write_h5mu(str(h5mu_path))
    return mdata, h5mu_path


@pytest.fixture
def synthetic_mudata(tmp_path):
    return _make_mudata(tmp_path)


def test_basic(run_component, tmp_path):
    """Proportions are computed and stored in uns and obsm."""
    mdata, h5mu_path = _make_mudata(tmp_path, n_donors=3, n_subpops=4, cells_per_group=10)
    output_path = tmp_path / "output.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
    ])

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    adata = result.mod["rna"]

    # uns["proportions"] must be a dict
    assert "proportions" in adata.uns
    prop_dict = adata.uns["proportions"]
    assert isinstance(prop_dict, dict)

    # Reconstruct DataFrame from dict and check shape
    # uns dict is column-first (subpop → donor → value); pd.DataFrame restores original shape
    prop_df = pd.DataFrame(prop_dict)
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
    mdata, h5mu_path = _make_mudata(tmp_path, n_donors=2, n_subpops=n_subpops, cells_per_group=5)
    output_path = tmp_path / "output_uniform.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
    ])

    result = read_h5mu(str(output_path))
    prop_dict = result.mod["rna"].uns["proportions"]
    prop_df = pd.DataFrame(prop_dict)
    expected = 1.0 / n_subpops
    np.testing.assert_allclose(prop_df.values, expected, atol=1e-10)


def test_custom_column_names(run_component, tmp_path):
    """Component works with non-default obs column names."""
    rng = np.random.default_rng(0)
    obs = pd.DataFrame({
        "sample": ["s1"] * 10 + ["s2"] * 10,
        "cell_class": (["A"] * 5 + ["B"] * 5) * 2,
    }, index=[f"c{i}" for i in range(20)])
    X = rng.integers(0, 50, size=(20, 5)).astype(float)
    var = pd.DataFrame(index=[f"g{g}" for g in range(5)])
    adata = AnnData(X=X, obs=obs, var=var)
    mdata = MuData({"rna": adata})
    h5mu_path = tmp_path / "custom_cols.h5mu"
    mdata.write_h5mu(str(h5mu_path))
    output_path = tmp_path / "output_custom.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
        "--obs_participant_id", "sample",
        "--obs_subpopulation", "cell_class",
        "--uns_output", "my_props",
        "--obsm_output", "my_props",
    ])

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

    with pytest.raises(Exception):
        run_component([
            "--input", str(h5mu_path),
            "--output", str(output_path),
            "--obs_subpopulation", "nonexistent_column",
        ])


def test_output_compression(run_component, tmp_path):
    """Component writes compressed output without error."""
    mdata, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_compressed.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
        "--output_compression", "gzip",
    ])

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    assert "proportions" in result.mod["rna"].uns


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
