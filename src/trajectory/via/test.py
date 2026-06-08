import sys
import pytest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/trajectory/via/via",
    "resources_dir": "src/trajectory/via/",
    "config": "src/trajectory/via/config.vsh.yaml",
}
## VIASH END


def _make_input(tmp_path, seed=0, n_cells=200, n_clusters=5):
    """Synthetic h5mu with a 2-D PHATE-like embedding and Leiden cluster labels."""
    rng = np.random.default_rng(seed)

    # Linear pseudotime embedded in 10-D (VIA needs >2 dims for graph construction)
    pt = np.linspace(0, 1, n_cells)
    embedding = np.column_stack(
        [pt]
        + [pt * rng.uniform(0.1, 0.9) + rng.normal(0, 0.05, n_cells) for _ in range(9)]
    )

    # Cluster labels (equal-size bins along pseudotime)
    cluster_ids = (pt * n_clusters).astype(int).clip(0, n_clusters - 1)
    cluster_labels = [str(c) for c in cluster_ids]

    obs = pd.DataFrame(
        {"leiden": cluster_labels},
        index=[f"c{i}" for i in range(n_cells)],
    )
    adata = ad.AnnData(
        X=rng.integers(0, 50, (n_cells, 3)).astype("float32"),
        obs=obs,
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    adata.obsm["X_phate"] = embedding.astype("float32")

    path = tmp_path / "input.h5mu"
    mu.MuData({"rna": adata}).write_h5mu(str(path))
    return path


def test_basic_pseudotime(run_component, tmp_path):
    """Basic run: pseudotime stored in obs, via_graph stored in uns."""
    h5mu = _make_input(tmp_path)
    out = tmp_path / "out.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--root_user",
            "0",
            "--output",
            str(out),
        ]
    )

    assert out.is_file()
    adata = mu.read_h5mu(str(out)).mod["rna"]
    assert "via_pseudotime" in adata.obs.columns
    assert "via_graph" in adata.uns

    pt = adata.obs["via_pseudotime"].values
    assert np.isfinite(pt).all()
    assert pt.min() >= 0.0


def test_pseudotime_range(run_component, tmp_path):
    """Pseudotime values span a meaningful range (not all identical)."""
    h5mu = _make_input(tmp_path, seed=1)
    out = tmp_path / "out2.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--root_user",
            "0",
            "--output",
            str(out),
        ]
    )

    pt = mu.read_h5mu(str(out)).mod["rna"].obs["via_pseudotime"].values
    assert pt.max() - pt.min() > 0.01, "Pseudotime has no spread"


def test_custom_output_keys(run_component, tmp_path):
    """Custom obs/uns output key names are respected."""
    h5mu = _make_input(tmp_path, seed=2)
    out = tmp_path / "out3.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--root_user",
            "0",
            "--obs_pseudotime",
            "my_pseudotime",
            "--uns_graph",
            "my_graph",
            "--output",
            str(out),
        ]
    )

    adata = mu.read_h5mu(str(out)).mod["rna"]
    assert "my_pseudotime" in adata.obs.columns
    assert "my_graph" in adata.uns


def test_root_cluster_label(run_component, tmp_path):
    """root_user as a cluster label string runs without error."""
    h5mu = _make_input(tmp_path, seed=3)
    out = tmp_path / "out4.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--root_user",
            "0",  # cluster label "0" exists
            "--output",
            str(out),
        ]
    )

    assert "via_pseudotime" in mu.read_h5mu(str(out)).mod["rna"].obs.columns


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
