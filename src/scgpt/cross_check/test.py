import sys
import pytest
import subprocess
from mudata import read_h5mu
from openpipelinetestutils.asserters import assert_shape_equal
import re
import pandas as pd

## VIASH START
meta = {
    'executable': './target/docker/scgpt/cross_check/cross_check',
    'resources_dir': './resources_test/scgpt/',
    'config': './src/scgpt/cross_check/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "test_resources/Kim2020_Lung.h5mu"

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
    
    assert {"str_batch", "batch_id"}.issubset(set(output_mudata.mod["rna"].obs.columns)), "Batch columns were not added."
    assert {"gene_name", "id_in_vocab"}.issubset(set(output_mudata.mod["rna"].var.columns)), "Gene columns were not added."    
    pd.testing.assert_series_equal(output_mudata.mod["rna"].var["id_in_vocab"],
                                   pd.Series(1, index=output_mudata.mod["rna"].var.index, name="id_in_vocab"))
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
    
    assert {"str_batch", "batch_id"}.issubset(set(output_mudata.mod["rna"].obs.columns)), "Batch columns were not added."
    assert {"gene_name", "id_in_vocab"}.issubset(set(output_mudata.mod["rna"].var.columns)), "Gene columns were not added."    
    pd.testing.assert_series_equal(output_mudata.mod["rna"].var["id_in_vocab"],
                                   pd.Series(1, index=output_mudata.mod["rna"].var.index, name="id_in_vocab"))
    assert_shape_equal(output_mudata, input_mudata)
 
    # TODO: figure out why this fails - the removed columns stay there in the asserter
    # output_mudata.var = output_mudata.var.drop(columns=["id_in_vocab"])
    # output_mudata.mod["rna"].obs = output_mudata.mod["rna"].obs.drop(columns=["str_batch", "batch_id"])
    # output_mudata.mod["rna"].var = output_mudata.mod["rna"].var.drop(columns=["id_in_vocab"])
    # assert_annotation_objects_equal(output_mudata, input_mudata)


def test_cross_check_custom_params(run_component, random_path):
    output_path = random_path(extension="h5ad")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--model_dir", "src/scgpt/model",
        "--gene_name_layer", "dummy_layer",
    ]