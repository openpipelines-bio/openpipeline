import pytest
import subprocess
from mudata import read_h5mu
import re
import sys
from glob import glob
from openpipelinetestutils.asserters import assert_annotation_objects_equal
import torch

## VIASH START
meta = {
    'executable': './target/docker/scgpt/annotation/annotation',
    'resources_dir': './resources_test/scgpt/',
    'config': './src/scgpt/annotation/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "test_resources/Kim2020_Lung_preprocessed.h5mu"
input_gene_ids = meta["resources_dir"] + "test_resources/Kim2020_Lung_gene_ids.pt"
input_values = meta["resources_dir"] + "test_resources/Kim2020_Lung_values.pt"
input_padding_mask = meta["resources_dir"] + "test_resources/Kim2020_Lung_padding_mask.pt"
model = meta["resources_dir"] + "source/best_model.pt"
model_config = meta["resources_dir"] + "source/args.json"
model_vocab = meta["resources_dir"] + "source/vocab.json" 

@pytest.fixture
def input_mudata_subset_cpu_run(write_mudata_to_file):
    mudata = read_h5mu(input_path)
    mudata.mod["rna"] = mudata.mod["rna"][:100]
    return write_mudata_to_file(mudata)

@pytest.fixture
def input_gene_ids_subset_cpu_run(random_path):
    output_path = random_path(extension="pt")
    torch.save(torch.load(input_gene_ids)[:100], output_path)
    return output_path

@pytest.fixture
def input_values_subset_cpu_run(random_path):
    output_path = random_path(extension="pt")
    torch.save(torch.load(input_values)[:100], output_path)
    return output_path

def test_annotation(run_component,
                    random_h5mu_path,
                    input_mudata_subset_cpu_run,
                    input_gene_ids_subset_cpu_run,
                    input_values_subset_cpu_run):
    output_path = random_h5mu_path()

    args = [
        "--input", input_mudata_subset_cpu_run,
        "--output",  output_path,
        "--modality", "rna",
        "--input_gene_ids", input_gene_ids_subset_cpu_run,
        "--input_values", input_values_subset_cpu_run,
        "--input_padding_mask", input_padding_mask,
        "--model", model,
        "--model_config", model_config,
        "--model_vocab", model_vocab,
        "--gene_name_layer", "gene_name",
        "--batch_id_layer", "batch_id",
        "--predicted_cell_type_id", "predicted_cell_type",
        "--pad_token", "<pad>",
        # "--mask_ratio", 0,
        # "--mask_value", -1,
        # "--pad_value", -2,
        # "--n_cls", 8,
        # "--n_input_bins", 51,
        # "--batch_size", 64,
        "--output_compression", "gzip"
    ]
    run_component(args)
    
    output_mudata = read_h5mu(output_path)
    assert "predicted_cell_type" in output_mudata.mod["rna"].obs.columns
    
if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
