import sys
import pytest
from pathlib import Path
import anndata
import mudata
import numpy as np

## VIASH START
meta = {
    'executable': './target/docker/labels_transfer/xgboost/xgboost',
    'resources_dir': './resources_test/'
}
## VIASH END

reference_file = f"{meta['resources_dir']}/annotation_test_data/TS_Blood_filtered.h5ad"
input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def test_args(tmp_path, request):
    obsm_features, obs_targets, output_uns_parameters = request.param

    tempfile_reference_file = tmp_path / "reference.h5ad"
    tempfile_input_file = tmp_path / "input.h5mu"

    # read reference
    reference_adata = anndata.read_h5ad(reference_file)

    # generate reference obs targets
    for i, target in enumerate(obs_targets):
        class_names = [str(idx) for idx in range(i + 2)]  # e.g. ["0", "1", "2"], the higher the level, the more the classes
        reference_adata.obs[target] = np.random.choice(class_names, size=reference_adata.n_obs)

    # read input query
    input_mudata = mudata.read_h5mu(input_file)
    input_rna_adata = input_mudata.mod["rna"]
    
    # generate features
    reference_adata.obsm[obsm_features] = np.random.normal(size=(reference_adata.n_obs, 30))
    input_rna_adata.obsm[obsm_features] = np.random.normal(size=(input_rna_adata.n_obs, 30))

    # write files
    reference_adata.write_h5ad(str(tempfile_reference_file))
    input_mudata.write_h5mu(str(tempfile_input_file))

    return tempfile_reference_file, reference_adata, tempfile_input_file, input_rna_adata, obsm_features, obs_targets, output_uns_parameters


@pytest.mark.parametrize("test_args", [("X_integrated_scvi", ["celltype"], None), ("X_int", ["ann_level_1", "ann_level_2", "ann_level_3"], "lab_tran")], indirect=True)
def test_label_transfer(run_component, test_args):
    tempfile_reference_file, _, tempfile_input_file, _, obsm_features, obs_targets, output_uns_parameters = test_args

    args = [
        "--input", str(tempfile_input_file),
        "--modality", "rna",
        "--input_obsm_features", obsm_features,
        "--reference", str(tempfile_reference_file),
        "--reference_obsm_features", obsm_features,
        "--reference_obs_targets", ";".join(obs_targets),
        "--output", "output.h5mu",
        "--model_output", "model_one_class",
        "--use_gpu", "false",
        "--max_depth", "6"
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
        assert output_data.mod["rna"].uns[exp_uns][target].get("method") == "XGBClassifier", f"Wrong method in parameters\noutput: {output_data.mod['rna'].uns}"
        assert output_data.mod["rna"].uns[exp_uns][target].get("max_depth") == 6, f"Wrong number of neighbors in parameters\noutput: {output_data.mod['rna'].uns}"


@pytest.mark.parametrize("test_args", [("X_int", ["ann_level_1", "ann_level_2", "ann_level_3"], "lab_tran")], indirect=True)
def test_retraining(run_component, test_args, tmp_path):
    output_model = tmp_path / "model_retraining"
    output_path = tmp_path / "output.h5mu"
    output2_path = tmp_path / "output2.h5mu"
    tempfile_reference_file, _, tempfile_input_file, _, obsm_features, obs_targets, output_uns_parameters = test_args

    # Train first 2 targets
    args = [
        "--modality", "rna",
        "--input_obsm_features", obsm_features,
        "--reference", str(tempfile_reference_file),
        "--reference_obsm_features", obsm_features,
        "--model_output", str(output_model)]
    
    if output_uns_parameters is not None:
        args.extend(["--output_uns_parameters", output_uns_parameters])
    
    args1 = args + [
        "--input", str(tempfile_input_file),
        "--output", str(output_path),
        "--reference_obs_targets", ",".join(obs_targets[:2]),
        "--max_depth", "6"]
    run_component(args1)

    assert output_path.is_file()

    # Add more targets
    # Now the code should use 2 previously trained models,
    # and train only the remaining targets
    args2 = args + [
        "--input", str(output_path),
        "--output", str(output2_path),
        "--reference_obs_targets", ",".join(obs_targets),
        "--max_depth", "4"]
    run_component(args2)

    assert output2_path.is_file()

    output_data = mudata.read_h5mu(output2_path)

    for target in obs_targets:
        assert f"{target}_pred" in output_data.mod["rna"].obs, f"Predictions are missing from output\noutput: {output_data.mod['rna'].obs}"
        assert f"{target}_uncertainty" in output_data.mod["rna"].obs, f"Uncertainties are missing from output\noutput: {output_data.mod['rna'].obs}"
        assert output_uns_parameters in output_data.mod["rna"].uns, f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
        assert target in output_data.mod["rna"].uns[output_uns_parameters], f"Parameters are missing from output\noutput: {output_data.mod['rna'].uns}"
        assert output_data.mod["rna"].uns[output_uns_parameters][target].get("method") == "XGBClassifier", f"Wrong method in parameters\noutput: {output_data.mod['rna'].uns}"
    
    assert output_data.mod["rna"].uns[output_uns_parameters][obs_targets[0]].get("max_depth") == 6, f"Wrong number of neighbors in parameters\noutput: {output_data.mod['rna'].uns}"
    assert output_data.mod["rna"].uns[output_uns_parameters][obs_targets[1]].get("max_depth") == 6, f"Wrong number of neighbors in parameters\noutput: {output_data.mod['rna'].uns}"
    assert output_data.mod["rna"].uns[output_uns_parameters][obs_targets[2]].get("max_depth") == 4, f"Wrong number of neighbors in parameters\noutput: {output_data.mod['rna'].uns}"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))