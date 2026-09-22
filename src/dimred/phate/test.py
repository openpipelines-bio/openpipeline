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


def _branching_manifold(n_obs, n_pcs, rng):
    """Y-shaped manifold: one trunk splitting into two branches.

    PHATE is a trajectory-preserving embedding, so a branching structure is the
    kind of input it is meant to resolve - a closed circle or a Gaussian blob
    would let a broken embedding pass unnoticed. Returns (coordinates, branch label).
    """
    n_trunk = n_obs // 2
    n_branch = (n_obs - n_trunk) // 2
    n_branch_b = n_obs - n_trunk - n_branch

    trunk_t = np.linspace(0.0, 1.0, n_trunk)
    trunk = np.column_stack([trunk_t, np.zeros(n_trunk)])

    branch_a_t = np.linspace(0.0, 1.0, n_branch)
    branch_a = np.column_stack([1.0 + branch_a_t, branch_a_t])

    branch_b_t = np.linspace(0.0, 1.0, n_branch_b)
    branch_b = np.column_stack([1.0 + branch_b_t, -branch_b_t])

    coords = np.vstack([trunk, branch_a, branch_b])
    branch = np.array(
        ["trunk"] * n_trunk + ["a"] * n_branch + ["b"] * n_branch_b, dtype=object
    )

    # Lift the 2-D manifold into n_pcs dimensions and add measurement noise
    loadings = rng.standard_normal((2, n_pcs))
    embedded = coords @ loadings + rng.standard_normal((n_obs, n_pcs)) * 0.05
    return embedded, branch


def _make_mudata(tmp_path, n_obs=120, n_pcs=10, seed=42):
    """Synthetic MuData with a branching manifold in the X_pca obsm slot."""
    rng = np.random.default_rng(seed)
    X_pca, branch = _branching_manifold(n_obs, n_pcs, rng)

    # No X: the component only reads .obsm
    obs = pd.DataFrame(
        {"branch": branch}, index=[f"cell_{i}" for i in range(n_obs)]
    )
    adata = AnnData(obs=obs)
    adata.obsm["X_pca"] = X_pca.astype(np.float32)

    mdata = MuData({"rna": adata})
    h5mu_path = tmp_path / "phate_input.h5mu"
    mdata.write_h5mu(str(h5mu_path))
    return mdata, h5mu_path


def test_basic(run_component, tmp_path):
    """PHATE runs and stores embedding with correct shape."""
    mdata, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
        ]
    )

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    adata = result.mod["rna"]

    assert "X_phate" in adata.obsm, "X_phate not found in .obsm"
    assert adata.obsm["X_phate"].shape == (mdata.mod["rna"].n_obs, 2), (
        f"Unexpected shape: {adata.obsm['X_phate'].shape}"
    )


def test_custom_obsm_output(run_component, tmp_path):
    """Custom --obsm_output key is used."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_custom.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obsm_output",
            "X_phate_custom",
        ]
    )

    result = read_h5mu(str(output_path))
    assert "X_phate_custom" in result.mod["rna"].obsm


def test_n_components(run_component, tmp_path):
    """--n_components controls embedding dimensionality."""
    mdata, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_3d.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--n_components",
            "3",
        ]
    )

    result = read_h5mu(str(output_path))
    assert result.mod["rna"].obsm["X_phate"].shape == (mdata.mod["rna"].n_obs, 3)


def test_custom_obsm_input(run_component, tmp_path):
    """Component works with a non-default --obsm_input key."""
    rng = np.random.default_rng(1)
    n_obs = 60
    obs = pd.DataFrame(index=[f"c{i}" for i in range(n_obs)])
    adata = AnnData(obs=obs)
    adata.obsm["proportions"] = rng.random((n_obs, 5)).astype(np.float32)
    mdata = MuData({"rna": adata})
    h5mu_path = tmp_path / "props_input.h5mu"
    mdata.write_h5mu(str(h5mu_path))
    output_path = tmp_path / "output_props.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--obsm_input",
            "proportions",
        ]
    )

    result = read_h5mu(str(output_path))
    assert "X_phate" in result.mod["rna"].obsm
    assert result.mod["rna"].obsm["X_phate"].shape == (n_obs, 2)


def test_fixed_t(run_component, tmp_path):
    """Fixed --t value runs without error."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_t.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--t",
            "10",
        ]
    )

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    assert "X_phate" in result.mod["rna"].obsm


def test_output_compression(run_component, tmp_path):
    """Compressed output is written correctly."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_gz.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--output_compression",
            "gzip",
        ]
    )

    assert output_path.is_file()
    result = read_h5mu(str(output_path))
    assert "X_phate" in result.mod["rna"].obsm


def test_missing_obsm_key_raises(run_component, tmp_path):
    """Missing --obsm_input raises an informative error."""
    _, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_err.h5mu"

    with pytest.raises(subprocess.CalledProcessError) as exc:
        run_component(
            [
                "--input",
                str(h5mu_path),
                "--output",
                str(output_path),
                "--obsm_input",
                "does_not_exist",
            ]
        )

    assert not output_path.is_file()
    assert "does_not_exist" in exc.value.stdout.decode("utf-8")


def test_branches_preserved(run_component, tmp_path):
    """The embedding must keep the two branches of the manifold apart.

    Without this the suite only checks the output shape, which a constant or
    shuffled embedding would also satisfy.
    """
    mdata, h5mu_path = _make_mudata(tmp_path)
    output_path = tmp_path / "output_branches.h5mu"

    run_component(
        [
            "--input",
            str(h5mu_path),
            "--output",
            str(output_path),
            "--n_components",
            "2",
        ]
    )

    result = read_h5mu(str(output_path))
    emb = result.mod["rna"].obsm["X_phate"]
    branch = result.mod["rna"].obs["branch"].to_numpy()

    tips_a = emb[branch == "a"][-10:]
    tips_b = emb[branch == "b"][-10:]

    within = 0.5 * (
        np.linalg.norm(tips_a - tips_a.mean(axis=0), axis=1).mean()
        + np.linalg.norm(tips_b - tips_b.mean(axis=0), axis=1).mean()
    )
    between = np.linalg.norm(tips_a.mean(axis=0) - tips_b.mean(axis=0))

    assert between > 2 * within, (
        f"Branch tips not separated in the PHATE embedding: "
        f"between={between:.4f}, within={within:.4f}"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
