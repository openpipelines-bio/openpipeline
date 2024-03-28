import pytest
import subprocess
from mudata import read_h5mu
import re
import sys

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

def test_annotation(run_component, random_path):
    output_path = random_path(extension="h5mu")
    print(f"\n\n\n\n\n\n\n\n{input_path}")
    print(meta)
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--input_gene_ids", input_gene_ids,
        "--input_values", input_values,
        "--input_padding_mask", input_padding_mask,
        "--model", model,
        "--model_config", model_config,
        "--model_vocab", model_vocab,
        # "--gene_name_layer", "gene_name",
        # "--batch_id_layer", "batch_id",
        # "--predicted_cell_type_id", "predicted_cell_type",
        # "--pad_token", "<pad>",
        # "--mask_ratio", 0,
        # "--mask_value", -1,
        # "--pad_value", -2,
        # "--n_cls", 8,
        # "--n_input_bins", 51,
        # "--batch_size", 64,
        # "--output_compression", None  
    ]
    run_component(args)
    
    output_mudata = read_h5mu(output_path)
    assert True == False
    
if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
