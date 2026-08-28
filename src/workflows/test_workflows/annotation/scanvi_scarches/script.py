from mudata import read_h5mu
import sys
import pytest
import pandas
import scvi

##VIASH START
par = {"input": "harmony_knn/output.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    # Test output data
    input_mudata = read_h5mu(par["input"])
    expected_obsm = ["X_integrated_scanvi", "X_scanvi_umap"]
    expected_obs = ["scanvi_pred", "scanvi_probabilities"]
    expected_obsp = [
        "scanvi_integration_distances",
        "scanvi_integration_connectivities",
    ]
    expected_uns = ["scanvi_integration_neighbors"]

    assert "rna" in list(input_mudata.mod.keys()), "Input should contain rna modality."
    assert all(key in list(input_mudata.mod["rna"].obsm) for key in expected_obsm), (
        f"Input mod['rna'] obsm columns should be: {expected_obsm}, found: {input_mudata.mod['rna'].obsm.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].obs) for key in expected_obs), (
        f"Input mod['rna'] obs columns should be: {expected_obs}, found: {input_mudata.mod['rna'].obs.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].obsp) for key in expected_obsp), (
        f"Input mod['rna'] obsp columns should be: {expected_obsp}, found: {input_mudata.mod['rna'].obsp.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].uns) for key in expected_uns), (
        f"Input mod['rna'] uns columns should be: {expected_uns}, found: {input_mudata.mod['rna'].uns.keys()}."
    )

    assert input_mudata.mod["rna"].obs["scanvi_pred"].dtype == "category"
    assert pandas.api.types.is_float_dtype(
        input_mudata.mod["rna"].obs["scanvi_probabilities"].dtype
    )

    # Test output model
    registry = scvi.model.base.BaseModelClass.load_registry(par["model"])
    model_args = registry["setup_args"]
    assert registry["model_name"] == "SCANVI", (
        f"Expected model to be SCANVI, found: {registry['model_name']}."
    )
    assert model_args["unlabeled_category"] == "Unknown", (
        f"Expected unlabeled_category to be 'Unkown', found: {model_args['unlabeled_category']}."
    )
    assert model_args["labels_key"] == par["obs_target"], (
        f"Expected labels_key to be {par['obs_target']}, found: {model_args['labels_key']}."
    )
    assert model_args["batch_key"] == par["obs_batch_label"], (
        f"Expected batch_key to be {par['obs_batch_label']}, found: {model_args['batch_key']}."
    )
    assert model_args["categorical_covariate_keys"] == par["obs_covariate"], (
        f"Expected categorical_covariate_keys to be {par['obs_covariate']}, found: {model_args['categorical_covariate_keys']}."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
