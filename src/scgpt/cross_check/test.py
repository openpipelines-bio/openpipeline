import sys
import pytest
import subprocess
from anndata import read_h5ad
from openpipelinetestutils.asserters import assert_annotation_objects_equal
import re

## VIASH START
meta = {
    'executable': './target/docker/scgpt/cross_check/cross_check',
    'resources_dir': './resources_test/',
    'config': './src/scgpt/cross_check/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "Kim2020_Lung.h5ad"

def test_cross_check(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--output_compression", "gzip"
        "load_model_vocab", "true",
        "ori_batch_layer_name", "sample",
        "batch_id_layer", "batch_id",
    ]
    
    run_component(args)
    
    assert output_path.is_file(), "No output was created."
    output_mudata = read_h5ad(output_path)
    input_mudata = read_h5ad(input_path)