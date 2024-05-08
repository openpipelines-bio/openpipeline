import pytest
from mudata import read_h5mu
import sys
import numpy as np
import subprocess
import re

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
                             tmp_path):

    output_annotation_file = tmp_path / "Kim2020_Lung_subset_annotated.h5mu"
    args = [
        "--input", input_path,
        "--output",  output_annotation_file,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--model", model,
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--obs_batch_label", "sample",
        "--obs_predicted_cell_type", "predicted_cell_type",
        "--dsbn", "True"
    ]
    run_component(args)

    output_mudata = read_h5mu(output_annotation_file)
    output_adata = output_mudata.mod["rna"]
    assert "predicted_cell_type" in output_adata.obs.keys(), "predicted_cell_type is not present in anndata obs keys"
    assert not all(np.isnan(output_adata.obs["predicted_cell_type"])), "predicted cell types contain nan values"
    assert output_adata.obs["predicted_cell_type"].dtype == "int64", "predicted cell types are not integers"

    # run withou dsbn
    output_annotation_file_without_dsbn = tmp_path / "Kim2020_Lung_subset_annotated_no_dsbn.h5mu"
    args = [
        "--input", input_path,
        "--output",  output_annotation_file_without_dsbn,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--model", model,
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--obs_batch_label", "sample",
        "--obs_predicted_cell_type", "predicted_cell_type",
        "--dsbn", "False"
    ]
    run_component(args)
    # Read output file
    output_mdata_no_dsbn = read_h5mu(output_annotation_file_without_dsbn)
    output_adata_no_dsbn = output_mdata_no_dsbn.mod["rna"]

    # Assert that embeddings without dsbn are different
    assert not (output_adata.obs["predicted_cell_type"] == output_adata_no_dsbn.obs["predicted_cell_type"]).all(), "Cell type predictions with and without dsbn are the same"


def test_annotation_dsbn_without_batch_labels(run_component, tmp_path):
    output_annotation_labels_without_dsbn = tmp_path / "Kim2020_Lung_subset_annotated_labels_without_dsbn.h5mu"

    args = [
        "--input", input_path,
        "--output",  output_annotation_labels_without_dsbn,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--model", model,
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--obs_predicted_cell_type", "predicted_cell_type",
        "--dsbn", "True",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: When dsbn is set to True, you are required to provide batch labels \(obs_batch_labels\)\.",
        err.value.stdout.decode('utf-8'))


def test_annotation_non_existing_keys(run_component, tmp_path):

    output_annotation_dummy_values = tmp_path / "Kim2020_Lung_subset_annotated_dummy_key.h5mu"

    # Test for non-existing tokenized values key
    args_2 = [
        "--input", input_path,
        "--output",  output_annotation_dummy_values,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "dummy_values_tokenized",
        "--model", model,
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--obs_predicted_cell_type", "predicted_cell_type",
        "--obs_batch_label", "sample",
        "--dsbn", "True",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_2)
    assert re.search(
        r'KeyError: "The parameter \'dummy_values_tokenized\' provided for \'--obsm_tokenized_values\' could not be found in adata.obsm"',
        err.value.stdout.decode('utf-8'))


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
