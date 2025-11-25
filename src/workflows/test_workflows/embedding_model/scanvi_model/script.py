from mudata import read_h5mu
import sys
import pytest
import pandas

##VIASH START
par = {"input": "harmony_knn/output.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["input"])
    expected_obsm = ["X_integrated_scanvi", "X_integrated_scvi"]
    expected_obs = ["scanvi_pred", "scanvi_probabilities"]

    assert "rna" in list(input_mudata.mod.keys()), "Input should contain rna modality."
    assert all(key in list(input_mudata.mod["rna"].obsm) for key in expected_obsm), (
        f"Input mod['rna'] obsm columns should be: {expected_obsm}, found: {input_mudata.mod['rna'].obsm.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].obs) for key in expected_obs), (
        f"Input mod['rna'] obs columns should be: {expected_obs}, found: {input_mudata.mod['rna'].obs.keys()}."
    )
    assert input_mudata.mod["rna"].obs["scanvi_pred"].dtype == "category"
    assert pandas.api.types.is_float_dtype(
        input_mudata.mod["rna"].obs["scanvi_probabilities"].dtype
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
