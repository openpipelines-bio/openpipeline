import re
import subprocess
import pytest
from pathlib import Path
import anndata as ad
import mudata as mu
import numpy as np
from scipy.sparse import csr_matrix

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


def test_label_transfer(run_component, random_h5mu_path):

    args = [
        "--input", input_file,
        "--modality", "rna",
        "--reference", reference_file,
        "--reference_obs_targets", "cell_type",
        "--output", random_h5mu_path(),
        "--n_neighbors", "5"
    ]

    run_component(args)

    assert Path("output.h5mu").is_file()

    output_data = mu.read_h5mu("output.h5mu")

    assert "cell_type_pred" in output_data.mod["rna"].obs, f"Predictions cell_type_pred is missing from output\noutput: {output_data.mod['rna'].obs}"
    assert "cell_type_probability" in output_data.mod["rna"].obs, f"Uncertainties cell_type_probability is missing from output\noutput: {output_data.mod['rna'].obs}"


@pytest.mark.parametrize("weights", ["uniform", "distance", "gaussian"])
def test_label_transfer_prediction_columns(run_component, weights, random_h5mu_path):

    output = random_h5mu_path()

    args = [
        "--input", input_file,
        "--modality", "rna",
        "--reference", reference_file,
        "--reference_obs_targets", "cell_type",
        "--weights", weights,
        "--output", output,
        "--output_obs_probability", "test_probability",
        "--output_obs_predictions", "test_prediction",
        "--n_neighbors", "5"
    ]

    run_component(args)

    assert Path(output).is_file()

    output_data = mu.read_h5mu(output)

    assert "test_prediction" in output_data.mod["rna"].obs, f"Predictions test_prediction is missing from output\noutput: {output_data.mod['rna'].obs}"
    assert "test_probability" in output_data.mod["rna"].obs, f"Uncertainties test_probability is missing from output\noutput: {output_data.mod['rna'].obs}"


def test_label_transfer_prediction_precomputed_neighbor_graph(run_component, random_h5mu_path):

    output = random_h5mu_path()

    # Add mock distance matrix to obsm slot
    reference_mdata = mu.read_h5mu(reference_file)
    ref_distances = np.random.rand(400, 400)
    ref_distances[ref_distances < 0.5] = 0
    ref_distances = csr_matrix(ref_distances)
    reference_mdata.mod["rna"].obsm["distances"] = ref_distances
    reference_mdata.write_h5mu(reference_file)

    query_mdata = mu.read_h5mu(input_file)
    query_distances = np.random.rand(713, 400)
    query_distances[query_distances < 0.5] = 0
    query_distances = csr_matrix(query_distances)
    query_mdata.mod["rna"].obsm["distances"] = query_distances
    query_mdata.write_h5mu(input_file)

    args = [
        "--input", input_file,
        "--modality", "rna",
        "--reference", reference_file,
        "--reference_obs_targets", "cell_type",
        "--output", output,
        "--input_obsm_distances", "distances",
        "--reference_obsm_distances", "distances",
        "--output_obs_probability", "test_probability",
        "--output_obs_predictions", "test_prediction",
        "--n_neighbors", "5"
    ]

    run_component(args)

    assert Path(output).is_file()

    output_data = mu.read_h5mu(output)

    assert "test_prediction" in output_data.mod["rna"].obs, f"Predictions test_prediction is missing from output\noutput: {output_data.mod['rna'].obs}"
    assert "test_probability" in output_data.mod["rna"].obs, f"Uncertainties test_probability is missing from output\noutput: {output_data.mod['rna'].obs}"


def test_raises_distance_matrix_dimensions(run_component, random_h5mu_path):

    output = random_h5mu_path()

    reference_mdata = mu.read_h5mu(reference_file)
    ref_distances = np.random.rand(400, 100)
    ref_distances[ref_distances < 0.5] = 0
    ref_distances = csr_matrix(ref_distances)
    reference_mdata.mod["rna"].obsm["distances"] = ref_distances
    reference_mdata.write_h5mu(reference_file)

    query_mdata = mu.read_h5mu(input_file)
    query_distances = np.random.rand(713, 400)
    query_distances[query_distances < 0.5] = 0
    query_distances = csr_matrix(query_distances)
    query_mdata.mod["rna"].obsm["distances"] = query_distances
    query_mdata.write_h5mu(input_file)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--modality", "rna",
            "--reference", reference_file,
            "--reference_obs_targets", "cell_type",
            "--output", output,
            "--input_obsm_distances", "distances",
            "--reference_obsm_distances", "distances",
            "--output_obs_probability", "test_probability",
            "--output_obs_predictions", "test_prediction",
            "--n_neighbors", "5"
        ])
    assert re.search(
        r"ValueError: The number of neighbors in the query and reference distance matrices do not match. Make sure both distance matrices contain distances to the reference dataset.",
        err.value.stdout.decode('utf-8')
        )


if __name__ == '__main__':
    exit(pytest.main([__file__]))
