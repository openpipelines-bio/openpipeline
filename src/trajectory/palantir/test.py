import sys
import pytest
import numpy as np
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/trajectory/palantir/palantir",
    "resources_dir": "src/trajectory/palantir/",
    "config": "src/trajectory/palantir/config.vsh.yaml",
}
## VIASH END


def _make_mudata(tmp_path, n_cells=300, n_dims=10, seed=42):
    """Synthetic MuData with a linear trajectory in X_pca_integrated.

    Cells are arranged along a 1-D pseudotime axis embedded in n_dims dimensions
    with small noise, giving Palantir a clear gradient to recover.
    """
    rng = np.random.default_rng(seed)
    pt = np.linspace(0, 1, n_cells)

    # Linear trajectory + noise across all dims
    embedding = np.column_stack(
        [pt]
        + [
            pt * rng.uniform(0.5, 1.0) + rng.normal(0, 0.05, n_cells)
            for _ in range(n_dims - 1)
        ]
    ).astype("float32")

    # Minimal gene matrix (Palantir doesn't use X, just the embedding)
    X = rng.negative_binomial(2, 0.5, (n_cells, 50)).astype("float32")

    obs = {
        "cluster": np.repeat(["A", "B", "C", "D", "E"], n_cells // 5),
    }
    obs["cluster"] = obs["cluster"][:n_cells]

    adata = ad.AnnData(
        X=X,
        obs={"cluster": obs["cluster"]},
    )
    adata.obs.index = [f"cell_{i:04d}" for i in range(n_cells)]
    adata.obsm["X_pca_integrated"] = embedding

    mdata = mu.MuData({"rna": adata})
    path = tmp_path / "palantir_input.h5mu"
    mdata.write_h5mu(str(path))
    return path


def test_palantir_start_cell_barcode(run_component, tmp_path):
    """Run Palantir using an explicit --start_cell barcode."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obsm_input",
            "X_pca_integrated",
            "--start_cell",
            "cell_0000",
            "--num_waypoints",
            "50",
            "--n_components",
            "5",
            "--knn",
            "10",
            "--output",
            str(output),
        ]
    )

    assert output.is_file(), "Output h5mu not created"

    adata = mu.read_h5mu(str(output)).mod["rna"]

    assert "palantir_pseudotime" in adata.obs.columns
    assert "palantir_entropy" in adata.obs.columns
    assert "palantir_fate_probabilities" in adata.obsm
    assert "palantir_waypoints" in adata.uns

    pt = adata.obs["palantir_pseudotime"]
    assert pt.notna().all(), "Pseudotime contains NaN values"
    assert (pt >= 0).all() and (pt <= 1).all(), "Pseudotime outside [0, 1]"
    assert len(adata.uns["palantir_waypoints"]) > 0, "Waypoints list is empty"


def test_palantir_start_cell_cluster(run_component, tmp_path):
    """Start cell is resolved automatically from a cluster label."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "output_cluster.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obsm_input",
            "X_pca_integrated",
            "--start_cell_cluster",
            "A",
            "--start_cell_obs_key",
            "cluster",
            "--num_waypoints",
            "50",
            "--n_components",
            "5",
            "--knn",
            "10",
            "--output",
            str(output),
        ]
    )

    assert output.is_file()
    adata = mu.read_h5mu(str(output)).mod["rna"]
    assert "palantir_pseudotime" in adata.obs.columns
    assert adata.obs["palantir_pseudotime"].notna().all()


def test_palantir_custom_output_keys(run_component, tmp_path):
    """Custom obs/obsm/uns output key names are respected."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "output_custom.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obsm_input",
            "X_pca_integrated",
            "--start_cell",
            "cell_0000",
            "--num_waypoints",
            "50",
            "--n_components",
            "5",
            "--knn",
            "10",
            "--obs_pseudotime",
            "my_pseudotime",
            "--obs_entropy",
            "my_entropy",
            "--obsm_fate_probabilities",
            "my_fate_probs",
            "--uns_waypoints",
            "my_waypoints",
            "--output",
            str(output),
        ]
    )

    adata = mu.read_h5mu(str(output)).mod["rna"]
    assert "my_pseudotime" in adata.obs.columns
    assert "my_entropy" in adata.obs.columns
    assert "my_fate_probs" in adata.obsm
    assert "my_waypoints" in adata.uns


def test_palantir_input_preserved(run_component, tmp_path):
    """Cell barcodes and gene names are unchanged after Palantir."""
    input_path = _make_mudata(tmp_path)
    output = tmp_path / "output_preserved.h5mu"

    run_component(
        [
            "--input",
            str(input_path),
            "--obsm_input",
            "X_pca_integrated",
            "--start_cell",
            "cell_0000",
            "--num_waypoints",
            "50",
            "--n_components",
            "5",
            "--knn",
            "10",
            "--output",
            str(output),
        ]
    )

    orig = mu.read_h5mu(str(input_path)).mod["rna"]
    out = mu.read_h5mu(str(output)).mod["rna"]

    np.testing.assert_array_equal(orig.obs_names, out.obs_names)
    np.testing.assert_array_equal(orig.var_names, out.var_names)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
