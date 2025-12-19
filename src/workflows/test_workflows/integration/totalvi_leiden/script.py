from mudata import read_h5mu
import sys
import pytest

##VIASH START
par = {"input": "harmony_knn/output.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["input"])
    
    expected_mod = ["rna", "prot"]
    expected_obsm_rna = ["X_integrated_totalvi", "X_totalvi_normalized_rna"]
    expected_obsp_rna = ["totalvi_integration_connectivities", "totalvi_integration_distances"]
    expected_obsm_prot = ["X_totalvi_normalized_protein"]
    

    assert all(key in list(input_mudata.mod) for key in expected_mod), (
        f"Input modalities should be: {expected_mod}, found: {input_mudata.mod.keys()}."
    )
    assert all(key in list(input_mudata.mod["rna"].obsm) for key in expected_obsm_rna), (
        f"Input mod['rna'] obsm columns should be: {expected_obsm_rna}, found: {input_mudata.mod['rna'].obsm.keys()}."
    )
    assert all(key in list(input_mudata.mod["prot"].obsm) for key in expected_obsm_prot), (
        f"Input mod['prot'] obsm columns should be: {expected_obsm_prot}, found: {input_mudata.mod['prot'].obsm.keys()}."
    )

if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
