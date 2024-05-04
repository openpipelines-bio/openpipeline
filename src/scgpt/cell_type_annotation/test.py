import pytest
from mudata import read_h5mu
import sys

## VIASH START
meta = {
    'executable': './target/docker/scgpt/annotation/annotation',
    'resources_dir': './resources_test/scgpt/',
    'config': './src/scgpt/annotation/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "Kim2020_Lung_subset_tokenized.h5mu"
model = meta["resources_dir"] + "best_model.pt"
model_config = meta["resources_dir"] + "args.json"
model_vocab = meta["resources_dir"] + "vocab.json" 

def test_cell_type_inference(run_component,
                             random_h5mu_path):
    
    output_path = random_h5mu_path()
    args = [
        "--input", input_path,
        "--output",  output_path,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--model", model,
        "--model_vocab", model_vocab,
        "--obs_batch_label", "sample",
        "--predicted_cell_type_id", "predictions",
        "--pad_token", "<pad>",
        "--dsbn", "True",
        "--pad_value", "-2",
        "--n_cls", "8",
        "--n_input_bins", "51",
        "--batch_size", "64"
    ]
    run_component(args)

    output_mudata = read_h5mu(output_path)
    assert "predictions" in output_mudata.mod["rna"].obs.columns
    assert output_mudata.mod["rna"].obs["predictions"].isna().sum() == 0
    assert output_mudata.mod["rna"].obs["predictions"].dtype == "int64"
    

if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))