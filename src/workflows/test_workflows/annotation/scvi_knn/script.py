from mudata import read_h5mu
import shutil
import os
import sys
from pathlib import Path
import pytest

##VIASH START
par = {"input": "scvi_knn/output.h5mu"}

meta = {"resources_dir": "resources_test"}
##VIASH END


def test_run():
    input_mudata = read_h5mu(par["input"])
    expected_obsm = ["X_integrated_scvi", "X_leiden_scvi_umap"]
    expected_obs = ["cell_type_pred", "cell_type_probability"]
    expected_obsp = ["scvi_integration_connectivities", "scvi_integration_distances"]

    assert "rna" in list(input_mudata.mod.keys()), "Input should contain rna modality."
    assert all(
        key in list(input_mudata.mod["rna"].obsm) for key in expected_obsm
    ), f"Input mod['rna'] obsm columns should be: {expected_obsm}, found: {input_mudata.mod['rna'].obsm.keys()}."
    assert all(
        key in list(input_mudata.mod["rna"].obs) for key in expected_obs
    ), f"Input mod['rna'] obs columns should be: {expected_obs}, found: {input_mudata.mod['rna'].obs.keys()}."
    assert all(
        key in list(input_mudata.mod["rna"].obsp) for key in expected_obsp
    ), f"Input mod['rna'] obsp columns should be: {expected_obsp}, found: {input_mudata.mod['rna'].obsp.keys()}."


if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    shutil.copyfile(
        os.path.join(meta["resources_dir"], "openpipelinetestutils", "conftest.py"),
        os.path.join(HERE_DIR, "conftest.py"),
    )
    sys.exit(pytest.main(["--import-mode=importlib"]))