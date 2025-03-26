import pytest
import sys
import mudata as mu
from scipy.sparse import issparse

## VIASH START
meta = {
    "resources_dir": "resources_test",
    "executable": "./target/docker/scgpt/binning/binning",
    "temp_dir": "tmp",
    "config": "./target/docker/scgpt/binning/.config.vsh.yaml",
}
## VIASH END


def test_binning(run_component, tmp_path):
    input_file_path = f"{meta['resources_dir']}/Kim2020_Lung_subset_preprocessed.h5mu"
    output_file_path = tmp_path / "Kim2020_Lung_subset_binned.h5mu"

    run_component(
        [
            "--input",
            input_file_path,
            "--modality",
            "rna",
            "--output_obsm_binned_counts",
            "binned_counts",
            "--n_input_bins",
            "51",
            "--var_input",
            "filter_with_hvg",
            "--output",
            output_file_path,
        ]
    )

    # Read output file
    output_mdata = mu.read(output_file_path)
    output_adata = output_mdata.mod["rna"]

    # Check presence of binning layers
    assert {"bin_edges", "binned_counts"}.issubset(
        output_adata.obsm.keys()
    ), "Binning obsm fields were not added."

    # Check bin edges
    bin_edges = output_adata.obsm["bin_edges"]
    assert all(bin_edges[:, 0] == 0)
    assert bin_edges.shape[1] == 51
    assert all(all(i >= 0) for i in bin_edges)

    # Check binned values
    binned_values = output_adata.obsm["binned_counts"]
    assert issparse(binned_values)
    assert (binned_values.data <= 51).all(axis=None)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
