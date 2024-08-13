import sys
import os
import pytest
import subprocess
import re
import mudata as mu
import scanpy as sc
import anndata as ad
from openpipelinetestutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "resources_dir": "resources_test"
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


def test_train_test(run_component, random_h5mu_path):
    output_train = random_h5mu_path()
    output_test = random_h5mu_path()
    
    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--test_size", "0.2",
        "--output_train", output_train,
        "--output_test", output_test,
    ])
    
    assert os.path.exists(output_train), "train file does not exist"
    assert os.path.exists(output_test), "test file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    train_mudata = mu.read_h5mu(output_train)
    test_mudata = mu.read_h5mu(output_test)
    
    assert list(train_mudata.mod.keys()) == list(test_mudata.mod.keys()) == ["rna"]
    
    assert train_mudata.mod["rna"].n_obs + test_mudata.mod["rna"].n_obs == input_mudata.mod["rna"].n_obs, \
        "train and test data do not sum up to input data"
        
    assert abs(train_mudata.mod["rna"].n_obs - round(input_mudata.mod["rna"].n_obs * 0.8)) <= 1, \
        "train data has wrong size"
    assert abs(test_mudata.mod["rna"].n_obs - round(input_mudata.mod["rna"].n_obs * 0.2)) <= 1, \
        "test data has wrong size"


def test_train_val_test(run_component, random_h5mu_path):
    output_train = random_h5mu_path()
    output_val = random_h5mu_path()
    output_test = random_h5mu_path()
    
    run_component([
        "--input", input_file,
        "--modality", "rna",
        "--test_size", "0.2",
        "--val_size", "0.1",
        "--output_train", output_train,
        "--output_val", output_val,
        "--output_test", output_test,
    ])
    
    assert os.path.exists(output_train), "train file does not exist"
    assert os.path.exists(output_val), "val file does not exist"
    assert os.path.exists(output_test), "test file does not exist"

    input_mudata = mu.read_h5mu(input_file)
    train_mudata = mu.read_h5mu(output_train)
    val_mudata = mu.read_h5mu(output_val)
    test_mudata = mu.read_h5mu(output_test)
    
    assert list(train_mudata.mod.keys()) == list(val_mudata.mod.keys()) == list(test_mudata.mod.keys()) == ["rna"]
    
    assert train_mudata.mod["rna"].n_obs + val_mudata.mod["rna"].n_obs + test_mudata.mod["rna"].n_obs == input_mudata.mod["rna"].n_obs, \
        "train, val and test data do not sum up to input data"
        
    assert abs(train_mudata.mod["rna"].n_obs - round(input_mudata.mod["rna"].n_obs * 0.7)) <= 1, \
        "train data has wrong size"
    assert abs(val_mudata.mod["rna"].n_obs - round(input_mudata.mod["rna"].n_obs * 0.1)) <= 1, \
        "val data has wrong size"
    assert abs(test_mudata.mod["rna"].n_obs - round(input_mudata.mod["rna"].n_obs * 0.2)) <= 1, \
        "test data has wrong size"

def test_raise_test_val_size(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--modality", "rna",
            "--test_size", "0.9",
            "--val_size", "0.5",
            "--output_train", "train.h5mu",
            "--output_val", "val.h5mu",
            "--output_test", "test.h5mu",
        ])
    
    assert re.search(r"Sum of test_size and val_size must not exceed 1.",
        err.value.stdout.decode('utf-8'))
    
    
def test_raise_invalid_val_out(run_component, random_h5mu_path):
    
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component([
            "--input", input_file,
            "--modality", "rna",
            "--test_size", "0.2",
            "--val_size", "0.1",
            "--output_train", "train.h5mu",
            "--output_test", "test.h5mu",
        ])
    
    assert re.search(r"Both --val_size and --output_val must be set to use validation set.",
        err.value.stdout.decode('utf-8'))

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))