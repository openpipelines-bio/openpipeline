import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest

## VIASH START
meta = {
    "executable": "./target/executable/filter/create_cell_masks/create_cell_masks",
    "resources_dir": "src/utils",
    "config": "./src/filter/create_cell_masks/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def input_h5mu():
    rng = np.random.default_rng(seed=1)
    n_obs = 100

    rna = ad.AnnData(
        X=rng.integers(0, 100, size=(n_obs, 10)),
        obs=pd.DataFrame(
            {
                "total_counts": rng.integers(0, 100, size=n_obs),
                "control_probe_counts": rng.integers(0, 20, size=n_obs),
            },
            index=[f"cell_{i}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(index=[f"gene_{i}" for i in range(10)]),
    )
    return mu.MuData({"rna": rna})


@pytest.fixture
def input_path(input_h5mu, random_h5mu_path):
    path = random_h5mu_path()
    input_h5mu.write(path)
    return path


def test_simple_execution(run_component, input_path, random_h5mu_path):
    output = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output,
            "--modality",
            "rna",
            "--filters",
            "total_counts:gt:10:rna,control_probe_counts:le:10:control,total_counts:le:50:rna",
        ]
    )
    assert output.is_file(), "output MuData file was not created"

    mdata = mu.read_h5mu(str(output))
    table = mdata["rna"]

    # Individual masks in obsm
    assert "cell_masks" in table.obsm, (
        "Individual masks not stored in obsm['cell_masks']"
    )
    masks_df = table.obsm["cell_masks"]
    assert isinstance(masks_df, pd.DataFrame), "Masks should be a DataFrame"
    assert len(masks_df.columns) == 6, (
        f"Expected 6 masks (3 individual, 3 combined), got {len(masks_df.columns)}"
    )
    assert "rna_total_counts_gt_10" in masks_df.columns
    assert "control_control_probe_counts_le_10" in masks_df.columns
    assert "rna_total_counts_le_50" in masks_df.columns
    assert "rna" in masks_df.columns
    assert "control" in masks_df.columns
    assert "overall" in masks_df.columns

    # Group masks in obs
    assert "cell_mask_rna" in table.obs.columns
    assert "cell_mask_control" in table.obs.columns
    assert "cell_mask" in table.obs.columns

    # Masks should be boolean
    assert table.obs["cell_mask_rna"].dtype == bool
    assert table.obs["cell_mask_control"].dtype == bool
    assert table.obs["cell_mask"].dtype == bool

    # Filters stored in uns
    assert "cell_filters" in table.uns
    filters_df = table.uns["cell_filters"]
    assert isinstance(filters_df, pd.DataFrame)
    assert len(filters_df) == 3

    filter_names = set(filters_df["name"].tolist())
    assert "rna_total_counts_gt_10" in filter_names
    assert "control_control_probe_counts_le_10" in filter_names
    assert "rna_total_counts_le_50" in filter_names

    filter_groups = set(filters_df["group"].tolist())
    assert "rna" in filter_groups
    assert "control" in filter_groups


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
