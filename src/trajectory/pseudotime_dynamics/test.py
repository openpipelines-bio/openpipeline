import sys
import pytest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/trajectory/pseudotime_dynamics/pseudotime_dynamics",
    "resources_dir": "src/trajectory/pseudotime_dynamics/",
    "config": "src/trajectory/pseudotime_dynamics/config.vsh.yaml",
}
## VIASH END


def _make_input(tmp_path, n_participants=15, cells_per=12, seed=0):
    """Synthetic h5mu with pseudotime + uns['proportions']."""
    rng = np.random.default_rng(seed)
    subpops = ["A", "B", "C", "D"]
    participant_ids = [f"donor_{i:02d}" for i in range(n_participants)]
    participant_pt = np.linspace(0.05, 0.95, n_participants)

    # Simple proportion trajectories
    def sig(x, x0):
        return 1 / (1 + np.exp(-10 * (x - x0)))

    props = np.column_stack(
        [
            sig(participant_pt, 0.3),
            1 - sig(participant_pt, 0.7),
            np.ones(n_participants) * 0.4,
            np.ones(n_participants) * 0.3,
        ]
    )
    props = np.clip(props + rng.normal(0, 0.02, props.shape), 0.01, None)
    props = (props.T / props.sum(axis=1)).T
    prop_df = pd.DataFrame(props, index=participant_ids, columns=subpops)

    obs_rows = []
    for i, pid in enumerate(participant_ids):
        for _ in range(cells_per):
            obs_rows.append(
                {
                    "participant_id": pid,
                    "subpopulation": rng.choice(subpops),
                    "palantir_pseudotime": float(
                        np.clip(participant_pt[i] + rng.normal(0, 0.02), 0, 1)
                    ),
                }
            )
    obs = pd.DataFrame(
        obs_rows, index=[f"c{i}" for i in range(n_participants * cells_per)]
    )
    adata = ad.AnnData(
        X=rng.integers(0, 100, (len(obs), 5)).astype("float32"),
        obs=obs,
        var=pd.DataFrame(index=[f"g{j}" for j in range(5)]),
    )
    adata.uns["proportions"] = prop_df.to_dict()

    path = tmp_path / "input.h5mu"
    mu.MuData({"rna": adata}).write_h5mu(str(path))
    return path, subpops


def test_basic(run_component, tmp_path):
    """Runs without error and stores uns['dynamics'] with expected keys."""
    h5mu, subpops = _make_input(tmp_path)
    out = tmp_path / "out.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--output",
            str(out),
            "--n_splines",
            "5",
        ]
    )

    assert out.is_file()
    mdata = mu.read_h5mu(str(out))
    dyn = mdata.mod["rna"].uns["dynamics"]

    assert set(dyn.keys()) == set(subpops), "Not all subpops in dynamics"
    for sp in subpops:
        assert "pseudotime_grid" in dyn[sp]
        assert "proportion_fitted" in dyn[sp]
        assert "peak_pseudotime" in dyn[sp]
        assert "r_squared" in dyn[sp]
        assert "p_value" in dyn[sp]
        assert len(dyn[sp]["pseudotime_grid"]) == 100  # default n_pseudotime_bins
        assert len(dyn[sp]["proportion_fitted"]) == 100


def test_custom_keys(run_component, tmp_path):
    """Custom obs/uns key names are respected."""
    h5mu, _ = _make_input(tmp_path, seed=1)
    # rename columns to use non-default names
    mdata = mu.read_h5mu(str(h5mu))
    adata = mdata.mod["rna"]
    adata.obs["my_pt"] = adata.obs["palantir_pseudotime"]
    adata.obs["my_pid"] = adata.obs["participant_id"]
    adata.uns["my_props"] = adata.uns["proportions"]
    path2 = tmp_path / "input2.h5mu"
    mdata.write_h5mu(str(path2))
    out = tmp_path / "out2.h5mu"

    run_component(
        [
            "--input",
            str(path2),
            "--obs_pseudotime",
            "my_pt",
            "--obs_participant_id",
            "my_pid",
            "--uns_proportions",
            "my_props",
            "--uns_output",
            "my_dynamics",
            "--n_splines",
            "5",
            "--output",
            str(out),
        ]
    )

    mdata_out = mu.read_h5mu(str(out))
    assert "my_dynamics" in mdata_out.mod["rna"].uns


def test_n_pseudotime_bins(run_component, tmp_path):
    """--n_pseudotime_bins controls the length of grid/fitted arrays."""
    h5mu, subpops = _make_input(tmp_path, seed=2)
    out = tmp_path / "out3.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--n_splines",
            "5",
            "--n_pseudotime_bins",
            "50",
            "--output",
            str(out),
        ]
    )

    dyn = mu.read_h5mu(str(out)).mod["rna"].uns["dynamics"]
    for sp in subpops:
        assert len(dyn[sp]["pseudotime_grid"]) == 50


def test_peak_pseudotime_range(run_component, tmp_path):
    """Peak pseudotime must be in [0, 1]."""
    h5mu, subpops = _make_input(tmp_path, seed=3)
    out = tmp_path / "out4.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--n_splines",
            "5",
            "--output",
            str(out),
        ]
    )

    dyn = mu.read_h5mu(str(out)).mod["rna"].uns["dynamics"]
    for sp in subpops:
        assert 0.0 <= dyn[sp]["peak_pseudotime"] <= 1.0, (
            f"peak_pseudotime out of range for {sp}"
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
