from mudata import read_h5mu
import sys
import pytest

##VIASH START
par = {"input": "harmony_knn/output.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["input"])
    expected_obs = ["celltypist_pred", "celltypist_probability"]
    expected_mod = par["expected_modalities"]

    assert all(key in list(input_mudata.mod) for key in expected_mod), (
        f"Input modalities should be: {expected_mod}, found: {input_mudata.mod.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].obs) for key in expected_obs), (
        f"Input mod['rna'] obs columns should be: {expected_obs}, found: {input_mudata.mod['rna'].obs.keys()}."
    )

    assert input_mudata.mod["rna"].obs["celltypist_pred"].dtype == "category", (
        "Cell type predictions should be of dtype category."
    )
    assert input_mudata.mod["rna"].obs["celltypist_probability"].dtype == "float64", (
        "Cell type probabilities should be of dtype float64."
    )

    if len(expected_mod) == 2:
        assert (
            input_mudata.mod[expected_mod[0]].shape[0]
            == input_mudata.mod[expected_mod[1]].shape[0]
        ), "Number of observations should be equal in all modalities."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
