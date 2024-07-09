import pytest
from mudata import read_h5mu
import sys
import torch
import numpy as np
import subprocess
import re


input_path = f'{meta["resources_dir"]}/Kim2020_Lung_subset_tokenized.h5mu'
model = f'{meta["resources_dir"]}/best_model.pt'
ft_model = f'{meta["resources_dir"]}/ft_best_model.pt'
model_config = f'{meta["resources_dir"]}/args.json'
model_vocab = f'{meta["resources_dir"]}/vocab.json'


def scgpt_to_ft_scgpt(scgpt_path, ft_scgpt_path, state_dict_key, mapper_key):
    f_model_dict = torch.load(scgpt_path, map_location="cpu")
    model_dict = {}
    model_dict[state_dict_key] = f_model_dict
    model_dict[mapper_key] = {k: str(k) for k in range(15)}
    torch.save(model_dict, ft_scgpt_path)


# Convert foundation model into fine-tuned model architecture:
# To be able to do a cell type label mapping, the model architecture needs to contain a class to label mapper dictionary
scgpt_to_ft_scgpt(model, ft_model, "model_state_dict", "id_to_class")


def test_cell_type_inference(run_component, tmp_path):

    output_annotation_file = tmp_path / "Kim2020_Lung_subset_annotated.h5mu"

    args = [
        "--input", input_path,
        "--output",  output_annotation_file,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--model", ft_model,
        "--finetuned_checkpoints_key", "model_state_dict",
        "--label_mapper_key", "id_to_class",
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--obs_batch_label", "sample",
        "--obs_predicted_cell_class", "predicted_cell_class",
        "--obs_predicted_cell_label", "predicted_cell_label",
        "--dsbn", "True"
    ]
    run_component(args)

    output_mudata = read_h5mu(output_annotation_file)
    output_adata = output_mudata.mod["rna"]
    assert "predicted_cell_class" in output_adata.obs.keys(), "predicted_cell_class is not present in anndata obs keys"
    assert not all(np.isnan(output_adata.obs["predicted_cell_class"])), "predicted cell classes contain nan values"
    assert output_adata.obs["predicted_cell_class"].dtype == "int64", "predicted cell classes are not integers"
    assert "predicted_cell_label" in output_adata.obs.keys(), "predicted_cell_type is not present in anndata obs keys"

    # run withou dsbn
    output_annotation_file_without_dsbn = tmp_path / "Kim2020_Lung_subset_annotated_no_dsbn.h5mu"
    args = [
        "--input", input_path,
        "--output",  output_annotation_file_without_dsbn,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--model", ft_model,
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--finetuned_checkpoints_key", "model_state_dict",
        "--label_mapper_key", "id_to_class",
        "--obs_batch_label", "sample",
        "--obs_predicted_cell_class", "predicted_cell_class",
        "--obs_predicted_cell_label", "predicted_cell_label",
        "--dsbn", "False"
    ]
    run_component(args)
    # Read output file
    output_mdata_no_dsbn = read_h5mu(output_annotation_file_without_dsbn)
    output_adata_no_dsbn = output_mdata_no_dsbn.mod["rna"]

    # Assert that embeddings without dsbn are different
    assert not (output_adata.obs["predicted_cell_class"] == output_adata_no_dsbn.obs["predicted_cell_class"]).all(), "Cell type predictions with and without dsbn are the same"


def test_annotation_dsbn_without_batch_labels(run_component, tmp_path):

    output_annotation_labels_without_dsbn = tmp_path / "Kim2020_Lung_subset_annotated_labels_without_dsbn.h5mu"

    args = [
        "--input", input_path,
        "--output",  output_annotation_labels_without_dsbn,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--model", ft_model,
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--finetuned_checkpoints_key", "model_state_dict",
        "--label_mapper_key", "id_to_class",
        "--obs_predicted_cell_class", "predicted_cell_class",
        "--obs_predicted_cell_label", "predicted_cell_label",
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
    args = [
        "--input", input_path,
        "--output",  output_annotation_dummy_values,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "dummy_values_tokenized",
        "--model", ft_model,
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--finetuned_checkpoints_key", "model_state_dict",
        "--label_mapper_key", "id_to_class",
        "--obs_predicted_cell_class", "predicted_cell_class",
        "--obs_predicted_cell_label", "predicted_cell_label",
        "--obs_batch_label", "sample",
        "--dsbn", "True",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r'KeyError: "The parameter \'dummy_values_tokenized\' provided for \'--obsm_tokenized_values\' could not be found in adata.obsm"',
        err.value.stdout.decode('utf-8'))

def test_checkpoint_architecture(run_component, tmp_path):

    output_dummy_model_key = tmp_path / "Kim2020_Lung_subset_annotated_dummy_key.h5mu"

    # Test for non-existing model file keys
    args = [
        "--input", input_path,
        "--output",  output_dummy_model_key,
        "--modality", "rna",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--model", ft_model,
        "--model_vocab", model_vocab,
        "--model_config", model_config,
        "--finetuned_checkpoints_key", "dummy_checkpoints_key",
        "--label_mapper_key", "id_to_class",
        "--obs_predicted_cell_class", "predicted_cell_class",
        "--obs_predicted_cell_label", "predicted_cell_label",
        "--obs_batch_label", "sample",
        "--dsbn", "True",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r'KeyError: "The key \'dummy_checkpoints_key\' provided for \'--finetuned_checkpoints_key\' could not be found in the provided --model file. The finetuned model file for cell type annotation requires valid keys for the checkpoints and the label mapper."',
        err.value.stdout.decode('utf-8'))


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
