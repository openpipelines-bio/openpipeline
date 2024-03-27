import pytest
import subprocess
from mudata import read_h5mu
import re

## VIASH START
meta = {
    'executable': './target/docker/scgpt/cross_check/cross_check',
    'resources_dir': './resources_test/scgpt/',
    'config': './src/scgpt/cross_check/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "test_resources/Kim2020_Lung_subset.h5mu"
vocab_path = meta["resources_dir"] + "source/vocab.json"

def test_cross_check(run_component, random_path):
    output_path = random_path(extension="h5mu")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--vocab_file", vocab_path,
    ]
    run_component(args)

    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(input_path)

    # Check added columns
    assert {"gene_name", "id_in_vocab"}.issubset(set(output_mudata.mod["rna"].var.columns)), "Gene columns were not added."    
    # Check if genes were filtered
    assert all(output_mudata.mod["rna"].var["id_in_vocab"] == 1), "Genes were not filtered."
    # Check if number of observations is the same
    assert output_mudata.mod["rna"].n_obs == input_mudata.mod["rna"].n_obs, "Number of observations changed."
    assert output_mudata.n_obs == input_mudata.n_obs, "Number of observations changed."

def test_cross_check_invalid_gene_layer_raises(run_component, random_path):
    output_path = random_path(extension="h5mu")
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--vocab_file", vocab_path,
        "--gene_name_layer", "dummy_gene",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(r"ValueError: gene name column 'dummy_gene' not found in .mod\['rna'\]\.obs\.",
                     err.value.stdout.decode('utf-8'))
    
    