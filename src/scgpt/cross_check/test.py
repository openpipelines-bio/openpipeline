import sys
import pytest
import subprocess
from anndata import read_h5ad
from openpipelinetestutils.asserters import assert_annotation_objects_equal
import re
import pandas as pd

## VIASH START
meta = {
    'executable': './target/docker/scgpt/cross_check/cross_check',
    'resources_dir': './resources_test/',
    'config': './src/scgpt/cross_check/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "Kim2020_Lung.h5ad"

def test_cross_check(run_component, random_path):
    output_path = random_path(extension="h5ad")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--model_dir", "src/scgpt/model",
        "--load_model_vocab", "true"
    ]
    
    run_component(args)
    
    assert output_path.is_file(), "No output was created."
    output_adata = read_h5ad(output_path)
    # input_adata = read_h5ad(input_path)
    
    assert {"str_batch", "batch_id"}.issubset(set(output_adata.obs.columns)), "Batch columns were not added."
    assert {"gene_name", "id_in_vocab"}.issubset(set(output_adata.var.columns)), "Gene columns were not added."    
    pd.testing.assert_series_equal(output_adata.var["id_in_vocab"],
                                   pd.Series(1, index=output_adata.var.index, name="id_in_vocab"))

    # TODO: add small synthetic data and compare expected to actual output

def test_cross_check_custom_params(run_component, random_path):
    output_path = random_path(extension="h5ad")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--model_dir", "src/scgpt/model",
        "--gene_name_layer", "dummy_layer",
    ]