import os
import sys
import pytest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/cluster/cellular_communities/cellular_communities",
    "resources_dir": "src/cluster/cellular_communities/",
    "config": "src/cluster/cellular_communities/config.vsh.yaml",
}
## VIASH END


def _make_input(tmp_path, seed=0, n_participants=16, cells_per=10):
    """Synthetic h5mu with uns['proportions'] and uns['dynamics']."""
    rng = np.random.default_rng(seed)
    subpops = ["A", "B", "C", "D", "E", "F"]
    n_cells = n_participants * cells_per
    pt = np.linspace(0.05, 0.95, n_participants)

    def sig(x, x0): return 1 / (1 + np.exp(-10 * (x - x0)))
    # 3 communities: (A,B), (C,D), (E,F)
    props = np.column_stack([sig(pt, 0.3), sig(pt, 0.32),
                              1 - sig(pt, 0.5), 1 - sig(pt, 0.48),
                              sig(pt, 0.7), sig(pt, 0.72)])
    props = np.clip(props + rng.normal(0, 0.01, props.shape), 0.01, None)
    props = (props.T / props.sum(axis=1)).T
    prop_df = pd.DataFrame(props, index=[f"d{i}" for i in range(n_participants)], columns=subpops)

    n_bins = 40
    grid = np.linspace(0, 1, n_bins)
    dynamics = {}
    for j, sp in enumerate(subpops):
        y = props[:, j][:n_bins]
        dynamics[sp] = {
            "pseudotime_grid": grid.tolist(),
            "proportion_fitted": y.tolist(),
            "peak_pseudotime": float(grid[np.argmax(y)]),
            "r_squared": 0.9,
            "p_value": 0.01,
        }

    obs_rows = [{"participant_id": f"d{i // cells_per}",
                 "subpopulation": rng.choice(subpops),
                 "palantir_pseudotime": float(np.clip(pt[i // cells_per] + rng.normal(0, 0.02), 0, 1))}
                for i in range(n_cells)]
    obs = pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(n_cells)])
    adata = ad.AnnData(X=rng.integers(0, 50, (n_cells, 3)).astype("float32"),
                       obs=obs,
                       var=pd.DataFrame(index=["g0", "g1", "g2"]))
    adata.uns["proportions"] = prop_df.to_dict()
    adata.uns["dynamics"] = dynamics

    path = tmp_path / "input.h5mu"
    mu.MuData({"rna": adata}).write_h5mu(str(path))
    return path, subpops


def test_basic_hierarchical(run_component, tmp_path):
    """Hierarchical clustering assigns all subpops to communities."""
    h5mu, subpops = _make_input(tmp_path)
    out = tmp_path / "out.h5mu"

    run_component([
        "--input", str(h5mu),
        "--n_communities", "3",
        "--output", str(out),
    ])

    assert out.is_file()
    mdata = mu.read_h5mu(str(out))
    adata = mdata.mod["rna"]

    assert "community_id" in adata.obs.columns
    assert "cellular_communities" in adata.uns

    comm_uns = adata.uns["cellular_communities"]
    assert set(comm_uns["subpopulation_communities"].keys()) == set(subpops)
    assert len(set(comm_uns["subpopulation_communities"].values())) == 3

    # every cell gets a community label
    assert adata.obs["community_id"].notna().all()


def test_n_communities_respected(run_component, tmp_path):
    """--n_communities controls the exact number of clusters produced."""
    h5mu, _ = _make_input(tmp_path, seed=1)
    out = tmp_path / "out2.h5mu"

    run_component([
        "--input", str(h5mu),
        "--n_communities", "2",
        "--output", str(out),
    ])

    comm = mu.read_h5mu(str(out)).mod["rna"].uns["cellular_communities"]
    assert len(set(comm["subpopulation_communities"].values())) == 2


def test_spectral_method(run_component, tmp_path):
    """Spectral clustering runs without error."""
    h5mu, _ = _make_input(tmp_path, seed=2)
    out = tmp_path / "out3.h5mu"

    run_component([
        "--input", str(h5mu),
        "--n_communities", "3",
        "--method", "spectral",
        "--output", str(out),
    ])

    assert out.is_file()
    assert "cellular_communities" in mu.read_h5mu(str(out)).mod["rna"].uns


def test_alpha_extremes(run_component, tmp_path):
    """alpha=1.0 (pure co-occurrence) and alpha=0.0 (pure dynamics) both run."""
    h5mu, _ = _make_input(tmp_path, seed=3)
    for alpha, suffix in [("1.0", "a"), ("0.0", "b")]:
        out = tmp_path / f"out_{suffix}.h5mu"
        run_component([
            "--input", str(h5mu),
            "--n_communities", "3",
            "--alpha", alpha,
            "--output", str(out),
        ])
        assert out.is_file()


def test_custom_output_keys(run_component, tmp_path):
    """Custom obs/uns output key names are respected."""
    h5mu, _ = _make_input(tmp_path, seed=4)
    out = tmp_path / "out_custom.h5mu"

    run_component([
        "--input", str(h5mu),
        "--n_communities", "3",
        "--obs_community_id", "my_community",
        "--uns_output", "my_communities",
        "--output", str(out),
    ])

    adata = mu.read_h5mu(str(out)).mod["rna"]
    assert "my_community" in adata.obs.columns
    assert "my_communities" in adata.uns


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
