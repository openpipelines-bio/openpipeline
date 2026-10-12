from mudata import read_h5mu
import numpy as np
import sys
import pytest

##VIASH START
par = {
    "input": "azimuth/output.h5mu",
    "modality": "rna",
    "expected_modalities": ["rna", "prot", "prot_imputed"],
    "expected_obs_predictions": ["azimuth_pred_celltype_l1"],
    "expected_obs_probability": ["azimuth_probability_celltype_l1"],
    "expected_obs_mapping_score": "azimuth_mapping_score",
    "expected_obsm": ["X_azimuth_umap", "X_integrated_azimuth"],
}

meta = {"resources_dir": "resources_test"}
##VIASH END


@pytest.fixture(scope="module")
def input_mudata():
    return read_h5mu(par["input"])


@pytest.fixture(scope="module")
def annotated(input_mudata):
    return input_mudata.mod[par["modality"]]


def test_modalities(input_mudata):
    assert sorted(input_mudata.mod) == sorted(par["expected_modalities"]), (
        f"Output modalities should be: {par['expected_modalities']}, found: {list(input_mudata.mod)}."
    )
    n_obs = input_mudata.mod[par["modality"]].n_obs
    for mod_name, mod in input_mudata.mod.items():
        assert mod.n_obs == n_obs, (
            f"Modality '{mod_name}' has {mod.n_obs} observations, expected {n_obs}."
        )


def test_predictions(annotated):
    expected_obs = (
        par["expected_obs_predictions"]
        + par["expected_obs_probability"]
        + [par["expected_obs_mapping_score"]]
    )
    missing = set(expected_obs) - set(annotated.obs.columns)
    assert not missing, (
        f"Missing .obs columns: {missing}, found: {list(annotated.obs.columns)}."
    )
    for col in par["expected_obs_predictions"]:
        assert annotated.obs[col].notna().all(), f"Column '{col}' contains NA values."
    for col in par["expected_obs_probability"] + [par["expected_obs_mapping_score"]]:
        values = annotated.obs[col]
        assert np.issubdtype(values.dtype, np.floating), (
            f"Column '{col}' should be a float, found: {values.dtype}."
        )
        # Prediction scores can exceed 1 by a floating point rounding error
        assert values.between(0, 1 + 1e-8).all(), (
            f"Column '{col}' should be between 0 and 1."
        )


def test_embeddings(annotated):
    for key in par["expected_obsm"]:
        assert key in annotated.obsm, (
            f"Missing .obsm key '{key}', found: {list(annotated.obsm)}."
        )
        assert annotated.obsm[key].shape[0] == annotated.n_obs


def test_counts_in_x(annotated):
    # The raw counts are stored in .X, and should not be duplicated in .layers
    assert annotated.X is not None, ".X should contain the raw counts."
    assert "X" not in annotated.layers, "Layer 'X' should have been removed."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
