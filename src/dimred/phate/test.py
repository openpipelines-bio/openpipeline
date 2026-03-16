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
    "executable": "./target/executable/dimred/phate/phate",
    "resources_dir": "./src/dimred/phate/",
    "cpus": 2,
    "config": "./src/dimred/phate/config.vsh.yaml",
}
## VIASH END


def _make_mudata(tmp_path, n_obs=80, n_pcs=10, seed=42):
    """Synthetic MuData with an X_pca obsm slot."""
    rng = np.random.default_rng(seed)
    # Simulate a simple manifold: points on a curve + noise
    t = np.linspace(0, 2 * np.pi, n_obs)
    coords = np.column_stack([np.sin(t), np.cos(t)])
    X_pca = coords @ rng.standard_normal((2, n_pcs)) + rng.standard_normal((n_obs, n_pcs)) * 0.1

    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=[f"gene_{g}" for g in range(5)])
    X = rng.integers(0, 100, size=(n_obs, 5)).astype(float)

    adata = AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_pca"] = X_pca.astype(np.float32)

    mdata = MuData({"rna": adata})
    h5mu_path = tmp_path / "phate_input.h5mu"
    mdata.write_h5mu(str(h5mu_path))
    return mdata, h5mu_path


def test_basic(run_component, tmp_path):
    """PHATE runs and stores embedding with correct shape."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
    ])

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    adata = result.mod["rna"]

    assert "X_phate" in adata.obsm, "X_phate not found in .obsm"
    assert adata.obsm["X_phate"].shape == (80, 2), \
        f"Unexpected shape: {adata.obsm['X_phate'].shape}"


def test_custom_obsm_output(run_component, tmp_path):
    """Custom --obsm_output key is used."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_custom.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
        "--obsm_output", "X_phate_custom",
    ])

    result = read_h5mu(str(output_path))
    assert "X_phate_custom" in result.mod["rna"].obsm


def test_n_components(run_component, tmp_path):
    """--n_components controls embedding dimensionality."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_3d.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
        "--n_components", "3",
    ])

    result = read_h5mu(str(output_path))
    assert result.mod["rna"].obsm["X_phate"].shape == (80, 3)


def test_custom_obsm_input(run_component, tmp_path):
    """Component works with a non-default --obsm_input key."""
    rng = np.random.default_rng(1)
    n_obs = 60
    obs = pd.DataFrame(index=[f"c{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=["g0"])
    adata = AnnData(
        X=rng.integers(0, 10, size=(n_obs, 1)).astype(float),
        obs=obs,
        var=var,
    )
    adata.obsm["proportions"] = rng.random((n_obs, 5)).astype(np.float32)
    mdata = MuData({"rna": adata})
    h5mu_path = tmp_path / "props_input.h5mu"
    mdata.write_h5mu(str(h5mu_path))
    output_path = tmp_path / "output_props.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
        "--obsm_input", "proportions",
    ])

    result = read_h5mu(str(output_path))
    assert "X_phate" in result.mod["rna"].obsm
    assert result.mod["rna"].obsm["X_phate"].shape == (n_obs, 2)


def test_fixed_t(run_component, tmp_path):
    """Fixed --t value runs without error."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_t.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
        "--t", "10",
    ])

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    assert "X_phate" in result.mod["rna"].obsm


def test_output_compression(run_component, tmp_path):
    """Compressed output is written correctly."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_gz.h5mu"

    run_component([
        "--input", str(h5mu_path),
        "--output", str(output_path),
        "--output_compression", "gzip",
    ])

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    assert "X_phate" in result.mod["rna"].obsm


def test_missing_obsm_key_raises(run_component, tmp_path):
    """Missing --obsm_input raises an informative error."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_err.h5mu"

    with pytest.raises(subprocess.CalledProcessError) as exc:
        run_component([
            "--input", str(h5mu_path),
            "--output", str(output_path),
            "--obsm_input", "does_not_exist",
        ])

    assert not output_path.is_file()
    assert "does_not_exist" in exc.value.stdout.decode("utf-8")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
