from mudata import read_h5mu
import numpy as np
import sys
import pytest

##VIASH START
par = {"input": "input.h5mu", "expected_modalities": ["rna", "prot"]}

meta = {"resources_dir": "resources_test"}
##VIASH END

EXPECTED_FEATURE_TYPES = {"rna": "Gene Expression", "prot": "Antibody Capture"}
VDJ_MODALITIES = ["vdj_t", "vdj_t_gd", "vdj_b"]


def test_run():
    input_mudata = read_h5mu(par["input"])
    expected_var = ["gene_name", "feature_type", "reference_file", "gene_ids"]
    expected_obs = ["run_id", "library_id", "cell_id"]

    assert set(input_mudata.mod.keys()) == set(par["expected_modalities"]), (
        f"Input should contain modalities {par['expected_modalities']}, "
        f"found: {list(input_mudata.mod.keys())}."
    )

    for modality, feature_type in EXPECTED_FEATURE_TYPES.items():
        if modality not in input_mudata.mod:
            continue
        adata = input_mudata.mod[modality]
        assert all(key in list(adata.var.columns) for key in expected_var), (
            f"Input mod['{modality}'] var columns should be: {expected_var}, found: {adata.var.keys()}."
        )
        assert all(key in list(adata.obs.columns) for key in expected_obs), (
            f"Input mod['{modality}'] obs columns should be: {expected_obs}, found: {adata.obs.keys()}."
        )
        assert np.array_equal(adata.var["feature_type"].unique(), [feature_type]), (
            f"Input mod['{modality}'] should only contain {feature_type} vars."
        )

    for modality in VDJ_MODALITIES:
        if modality not in input_mudata.mod:
            continue
        adata = input_mudata.mod[modality]
        assert "airr" in adata.obsm, (
            f"Input mod['{modality}'] should contain .obsm['airr']."
        )
        assert all(
            key in list(adata.obs.columns) for key in expected_obs + ["putative_cell"]
        ), (
            f"Input mod['{modality}'] obs columns should contain: {expected_obs + ['putative_cell']}, found: {adata.obs.keys()}."
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
