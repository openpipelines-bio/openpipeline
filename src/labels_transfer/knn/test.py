import pytest
from pathlib import Path
import anndata
import mudata
import numpy as np

## VIASH START
meta = {
    'executable': './target/executable/labels_transfer/knn/knn',
    'resources_dir': './resources_test/'
}
## VIASH END

reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"
input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"

@pytest.fixture
def test_args(tmp_path, request):
    obsm_features, obs_targets, output_uns_parameters = request.param

    # read reference
    tempfile_reference_file = tmp_path / "reference.h5ad"
    reference_adata = anndata.read_h5ad(reference_file)

    # generate reference obs targets
    for i, target in enumerate(obs_targets):
        class_names = [str(idx) for idx in range(i + 2)]  # e.g. ["0", "1", "2"], the higher the level, the more the classes
        reference_adata.obs[target] = np.random.choice(class_names, size=reference_adata.n_obs)

    # read input query
    tempfile_input_file = tmp_path / "input.h5mu"
    input_data = mudata.read_h5mu(input_file)
    adata = input_data.mod["rna"]
    
    # generate features
    reference_adata.obsm[obsm_features] = np.random.normal(size=(reference_adata.n_obs, 30))
    adata.obsm[obsm_features] = np.random.normal(size=(adata.n_obs, 30))

    # write files
    reference_adata.write(tempfile_reference_file.name)
    input_data.write(tempfile_input_file.name)

    return tempfile_reference_file, reference_adata, tempfile_input_file, adata, obsm_features, obs_targets, output_uns_parameters

@pytest.mark.parametrize("test_args", [("X_integrated_scvi", ["celltype"], None), ("X_int", ["ann_level_1", "ann_level_2", "ann_level_3"], "lab_tran")], indirect=True)
def test_label_transfer(run_component, test_args):
    tempfile_reference_file, reference_adata, tempfile_input_file, query_adata, obsm_features, obs_targets, output_uns_parameters = test_args

    args = [
        "--input", tempfile_input_file.name,
        "--modality", "rna",
        "--input_obsm_features", obsm_features,
        "--reference", tempfile_reference_file.name,
        "--reference_obsm_features", obsm_features,
        "--reference_obs_targets", ";".join(obs_targets),
        "--output", "output.h5mu",
        "--n_neighbors", "5"
    ]

    if output_uns_parameters is not None:
        args.extend(["--output_uns_parameters", output_uns_parameters])

    run_component(args)

    assert Path("output.h5mu").is_file()

    output_data = mudata.read_h5mu("output.h5mu")

    exp_uns = "labels_transfer" if output_uns_parameters is None else output_uns_parameters

    for target in obs_targets:
        assert f"{target}_pred" in output_data.mod["rna"].obs, f"Predictions are missing from output\noutput: {output_data.mod['rna'].obs}"
        assert f"{target}_uncertainty" in output_data.mod["rna"].obs, f"Uncertainties are missing from output\noutput: {output_data.mod['rna'].obs}"
        assert exp_uns in output_data.mod["rna"].uns, f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
        assert target in output_data.mod["rna"].uns[exp_uns], f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
        assert output_data.mod["rna"].uns[exp_uns][target].get("method") == "KNN_pynndescent", f"Wrong method in parameters\noutput: {output_data.mod['rna'].uns}"
        assert output_data.mod["rna"].uns[exp_uns][target].get("n_neighbors") == 5, f"Wrong number of neighbors in parameters\noutput: {output_data.mod['rna'].uns}"

if __name__ == '__main__':
    exit(pytest.main([__file__]))