import pytest
from pathlib import Path
import anndata
import mudata
import numpy as np

## VIASH START
meta = {
    'executable': './target/docker/labels_transfer/knn/knn',
    'resources_dir': './resources_test/'
}
## VIASH END

reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"
input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def test_label_transfer(run_component):

    args = [
        "--input", input_file,
        "--modality", "rna",
        "--reference", reference_file,
        "--reference_obs_targets", "cell_type"
        "--output", "output.h5mu",
        "--n_neighbors", "5"
    ]

    run_component(args)

    assert Path("output.h5mu").is_file()

    output_data = mudata.read_h5mu("output.h5mu")

    assert "cell_type_pred" in output_data.mod["rna"].obs, f"Predictions cell_type_pred is missing from output\noutput: {output_data.mod['rna'].obs}"
    assert "cell_type_probability" in output_data.mod["rna"].obs, f"Uncertainties cell_type_probability is missing from output\noutput: {output_data.mod['rna'].obs}"


if __name__ == '__main__':
    exit(pytest.main([__file__]))