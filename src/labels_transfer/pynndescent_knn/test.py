import pytest
from pathlib import Path
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    'resources_dir': './resources_test/'
}
## VIASH END

reference_h5ad_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"
# convert reference to h5mu
reference_adata = ad.read_h5ad(reference_h5ad_file)
reference_mdata = mu.MuData({"rna": reference_adata})
reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5mu"
reference_mdata.write_h5mu(reference_file)
input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def test_label_transfer(run_component):

    args = [
        "--input", input_file,
        "--modality", "rna",
        "--reference", reference_file,
        "--reference_obs_targets", "cell_type",
        "--output", "output.h5mu",
        "--n_neighbors", "5"
    ]

    run_component(args)

    assert Path("output.h5mu").is_file()

    output_data = mu.read_h5mu("output.h5mu")

    assert "cell_type_pred" in output_data.mod["rna"].obs, f"Predictions cell_type_pred is missing from output\noutput: {output_data.mod['rna'].obs}"
    assert "cell_type_probability" in output_data.mod["rna"].obs, f"Uncertainties cell_type_probability is missing from output\noutput: {output_data.mod['rna'].obs}"


def test_label_transfer_prediction_columns(run_component):

    args = [
        "--input", input_file,
        "--modality", "rna",
        "--reference", reference_file,
        "--reference_obs_targets", "cell_type",
        "--output", "output.h5mu",
        "--output_obs_probability", "test_probability",
        "--output_obs_predictions", "test_prediction",
        "--n_neighbors", "5"
    ]

    run_component(args)

    assert Path("output.h5mu").is_file()

    output_data = mu.read_h5mu("output.h5mu")

    assert "test_prediction" in output_data.mod["rna"].obs, f"Predictions test_prediction is missing from output\noutput: {output_data.mod['rna'].obs}"
    assert "test_probability" in output_data.mod["rna"].obs, f"Uncertainties test_probability is missing from output\noutput: {output_data.mod['rna'].obs}"


if __name__ == '__main__':
    exit(pytest.main([__file__]))
