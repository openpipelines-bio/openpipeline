import pytest
import subprocess
from mudata import read_h5mu
import sys
from openpipelinetestutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    'executable': './target/docker/scgpt/annotation/annotation',
    'resources_dir': './resources_test/scgpt/',
    'config': './src/scgpt/annotation/config.vsh.yaml'
}
## VIASH END

input_path = meta["resources_dir"] + "Kim2020_Lung_tokenized.h5mu"
model = meta["resources_dir"] + "best_model.pt"
model_config = meta["resources_dir"] + "args.json"
model_vocab = meta["resources_dir"] + "vocab.json" 

@pytest.fixture
def input_mudata_subset_cpu_run(write_mudata_to_file):
    mudata = read_h5mu(input_path)
    mudata.mod["rna"] = mudata.mod["rna"][:100]
    return write_mudata_to_file(mudata)

def test_cell_type_inference(run_component,
                             random_h5mu_path,
                             input_mudata_subset_cpu_run):
    output_path = random_h5mu_path()

    args = [
        "--input", input_mudata_subset_cpu_run,
        "--output",  output_path,
        "--modality", "rna",
        "--input_obsm_gene_tokens", "gene_id_tokens",
        "--input_obsm_tokenized_values", "values_tokenized",
        "--model", model,
        "--model_config", model_config,
        "--model_vocab", model_vocab,
        "--gene_name_layer", "gene_name",
        "--input_obs_batch_label", "sample",
        "--predicted_cell_type_id", "predictions",
        "--pad_token", "<pad>",
        "--DSBN", "True",
        "--pad_value", "-2",
        "--n_cls", "8",
        "--n_input_bins", "51",
        "--batch_size", "64"
    ]
    run_component(args)
    
    output_mudata = read_h5mu(output_path)
    input_mudata = read_h5mu(input_mudata_subset_cpu_run)
    assert "predictions" in output_mudata.mod["rna"].obs.columns
    assert output_mudata.mod["rna"].obs["predictions"].isna().sum() == 0
    assert output_mudata.mod["rna"].obs["predictions"].dtype == "int64"
    
    # print(output_mudata)
    # print(input_mudata)
    # output_mudata.mod["rna"].obs.drop("predictions",
    #                                   axis=1,
    #                                   inplace=True)
    
    # # print(input_mudata.mod["rna"].obs.columns)
    # # print(output_mudata.mod["rna"].obs.columns)
    # assert_annotation_objects_equal(output_mudata, input_mudata)
    
if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
