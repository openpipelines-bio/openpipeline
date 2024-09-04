from mudata import read_h5mu
import numpy as np
import shutil
import os
import sys
from pathlib import Path
import pytest

##VIASH START
par = {
    "input": "input.h5mu"
}

meta = {
    "resources_dir": "resources_test"
}
##VIASH END

def test_run():
    input_mudata = read_h5mu(par["input"])
    expected_layers = ["binned"]
    expected_obsm = ["gene_id_tokens", "values_tokenized", "padding_mask"]

    assert "rna" in list(input_mudata.mod.keys()), "Input should contain rna modality."

    assert all(key in list(input_mudata.mod["rna"].layers) for key in expected_layers), f"Input mod['rna'] var columns should be: {expected_layers}, found: {input_mudata.mod["rna"].var.keys()}."
    assert all(key in list(input_mudata.mod["rna"].obsm) for key in expected_obsm), f"Input mod['rna'] obs columns should be: {expected_obsm}, found: {input_mudata.mod["rna"].obs.keys()}."

    #TODO: check size based on hvg

if __name__ == "__main__":
    HERE_DIR = Path(__file__).resolve().parent
    shutil.copyfile(os.path.join(meta['resources_dir'], "openpipelinetestutils", "conftest.py"),
                    os.path.join(HERE_DIR, "conftest.py"))
    sys.exit(pytest.main(["--import-mode=importlib"]))