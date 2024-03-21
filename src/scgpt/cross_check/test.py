import pytest
import subprocess
from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_shape_equal
import re
import pandas as pd
import os

## VIASH START
meta = {
    'executable': './target/docker/scgpt/cross_check/cross_check',
    'resources_dir': './resources_test/scgpt/',
    'config': './src/scgpt/cross_check/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "test_resources/Kim2020_Lung.h5mu"

@pytest.fixture
def model_dir_missing_files(random_path):
    model_dir = random_path()
    os.mkdir(model_dir)
    args_json_path = os.path.join(model_dir, "args.json")
    vocab_json_path = os.path.join(model_dir, "vocab.json")
    open(args_json_path, 'w').close()
    open(vocab_json_path, 'w').close()
    return model_dir

def test_cross_check(run_component, random_path):
    output_path = random_path(extension="h5ad")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--model_dir", "resources_test/scgpt/source",
        "--load_model_vocab", "true"
    ]
    run_component(args)
    
    assert output_path.is_file(), "No output was created."
    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(input_path)
    
    # Check added columns
    assert {"str_batch", "batch_id"}.issubset(set(output_mudata.mod["rna"].obs.columns)), "Batch columns were not added."
    assert {"gene_name", "id_in_vocab"}.issubset(set(output_mudata.mod["rna"].var.columns)), "Gene columns were not added."    
    # Check if genes were filtered
    pd.testing.assert_series_equal(output_mudata.mod["rna"].var["id_in_vocab"],
                                   pd.Series(1, index=output_mudata.mod["rna"].var.index, name="id_in_vocab"))
    # Check if number of observations is the same
    assert output_mudata.mod["rna"].n_obs == input_mudata.mod["rna"].n_obs, "Number of observations changed."
    assert output_mudata.n_obs == input_mudata.n_obs, "Number of observations changed."

    
def test_cross_check_no_vocab(run_component, random_path):
    output_path = random_path(extension="h5ad")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--model_dir", "resources_test/scgpt/source",
        "--load_model_vocab", "false"
    ]
    run_component(args)
    
    assert output_path.is_file(), "No output was created."
    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(input_path)
    
    # Check added columns
    assert {"str_batch", "batch_id"}.issubset(set(output_mudata.mod["rna"].obs.columns)), "Batch columns were not added."
    assert {"gene_name", "id_in_vocab"}.issubset(set(output_mudata.mod["rna"].var.columns)), "Gene columns were not added."    
    # Check if genes were filtered
    pd.testing.assert_series_equal(output_mudata.mod["rna"].var["id_in_vocab"],
                                   pd.Series(1, index=output_mudata.mod["rna"].var.index, name="id_in_vocab"))
    
    assert_shape_equal(output_mudata, input_mudata)


def test_cross_check_custom_params(run_component, random_path):
    output_path = random_path(extension="h5ad")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--model_dir", "resources_test/scgpt/source",
        "--gene_name_layer", "dummy_layer",
        "--batch_id_layer", "dummy_batch_id"
    ]
    run_component(args)
    
    assert output_path.is_file(), "No output was created."
    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(input_path)
    
    # Check added columns
    assert {"dummy_batch_id"}.issubset(set(output_mudata.mod["rna"].obs.columns)), "Batch columns were not added."
    assert {"dummy_layer", "id_in_vocab"}.issubset(set(output_mudata.mod["rna"].var.columns)), "Gene columns were not added."
    # Check if genes were filtered
    pd.testing.assert_series_equal(output_mudata.mod["rna"].var["id_in_vocab"],
                                   pd.Series(1, index=output_mudata.mod["rna"].var.index, name="id_in_vocab"))
    # Check if number of observations is the same
    assert output_mudata.mod["rna"].n_obs == input_mudata.mod["rna"].n_obs, "Number of observations changed."
    assert output_mudata.n_obs == input_mudata.n_obs, "Number of observations changed."
    
    
def test_cross_check_invalid_batch_layer_raises(run_component, random_path):
    output_path = random_path(extension="h5ad")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--model_dir", "resources_test/scgpt/source",
        "--ori_batch_layer_name", "dummy_batch",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(r"ValueError: Batch column 'dummy_batch' not found in .mod\['rna'\]\.obs\.",
                     err.value.stdout.decode('utf-8'))

    
def test_cross_check_missing_model_files_raises(run_component, random_path, model_dir_missing_files):
    output_path = random_path(extension="h5ad")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--model_dir", model_dir_missing_files,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(f"FileNotFoundError: Model directory '/viash_automount{model_dir_missing_files}/' is missing the following required files: best_model.pt.",
                     err.value.stdout.decode('utf-8'))