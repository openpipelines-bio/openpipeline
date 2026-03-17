import os
import sys
import pytest
import numpy as np
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/trajectory/palantir/palantir",
    "resources_dir": "src/trajectory/palantir/",
    "config": "src/trajectory/palantir/config.vsh.yaml",
}
## VIASH END


@pytest.fixture(scope="module")
def atlas_data():
    """Load the test atlas once and return key information.
    Viash copies test_resources files directly into meta['resources_dir']."""
    atlas_path = os.path.join(meta["resources_dir"], "atlas.h5mu")
    mdata = mu.read_h5mu(atlas_path)
    adata = mdata.mod["rna"]
    # Pick a start cell from the first available subpopulation cluster
    first_cluster = adata.obs["subpopulation"].iloc[0]
    return {
        "path": atlas_path,
        "first_cluster": first_cluster,
        "first_cell": adata.obs_names[0],
        "n_cells": adata.n_obs,
    }


def test_palantir_start_cell_barcode(run_component, tmp_path, atlas_data):
    """Run Palantir using an explicit --start_cell barcode."""
    output = tmp_path / "output.h5mu"

    run_component(
        [
            "--input", atlas_data["path"],
            "--obsm_input", "X_pca_integrated",
            "--start_cell", atlas_data["first_cell"],
            "--num_waypoints", "50",
            "--n_components", "5",
            "--knn", "10",
            "--output", str(output),
        ]
    )

    assert output.is_file(), "Output h5mu not created"

    mdata = mu.read_h5mu(str(output))
    adata = mdata.mod["rna"]

    assert "palantir_pseudotime" in adata.obs.columns, \
        "palantir_pseudotime not in obs"
    assert "palantir_entropy" in adata.obs.columns, \
        "palantir_entropy not in obs"
    assert "palantir_fate_probabilities" in adata.obsm, \
        "palantir_fate_probabilities not in obsm"
    assert "palantir_waypoints" in adata.uns, \
        "palantir_waypoints not in uns"

    pt = adata.obs["palantir_pseudotime"]
    assert pt.notna().all(), "Pseudotime contains NaN values"
    assert (pt >= 0).all() and (pt <= 1).all(), \
        "Pseudotime values outside [0, 1]"

    assert len(adata.uns["palantir_waypoints"]) > 0, "Waypoints list is empty"


def test_palantir_custom_output_keys(run_component, tmp_path, atlas_data):
    """Custom obs/obsm/uns output key names are respected."""
    output = tmp_path / "output_custom.h5mu"

    run_component(
        [
            "--input", atlas_data["path"],
            "--start_cell", atlas_data["first_cell"],
            "--num_waypoints", "50",
            "--n_components", "5",
            "--knn", "10",
            "--obs_pseudotime", "my_pseudotime",
            "--obs_entropy", "my_entropy",
            "--obsm_fate_probabilities", "my_fate_probs",
            "--uns_waypoints", "my_waypoints",
            "--output", str(output),
        ]
    )

    mdata = mu.read_h5mu(str(output))
    adata = mdata.mod["rna"]

    assert "my_pseudotime" in adata.obs.columns
    assert "my_entropy" in adata.obs.columns
    assert "my_fate_probs" in adata.obsm
    assert "my_waypoints" in adata.uns


def test_palantir_input_preserved(run_component, tmp_path, atlas_data):
    """Cell barcodes and gene names are unchanged after Palantir."""
    output = tmp_path / "output_preserved.h5mu"

    run_component(
        [
            "--input", atlas_data["path"],
            "--start_cell", atlas_data["first_cell"],
            "--num_waypoints", "50",
            "--n_components", "5",
            "--knn", "10",
            "--output", str(output),
        ]
    )

    orig = mu.read_h5mu(atlas_data["path"])
    out = mu.read_h5mu(str(output))

    np.testing.assert_array_equal(
        orig.mod["rna"].obs_names, out.mod["rna"].obs_names
    )
    np.testing.assert_array_equal(
        orig.mod["rna"].var_names, out.mod["rna"].var_names
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
