import pytest
import subprocess
import re
import sys
import mudata as mu
import numpy as np


## VIASH START
meta = {
    "resources_dir": "resources_test",
}
## VIASH END

input = f"{meta['resources_dir']}/Kim2020_Lung_subset_tokenized.h5mu"
model_file = f"{meta['resources_dir']}/source/best_model.pt"
ft_model_file = f'{meta["resources_dir"]}/finetuned_model/best_model.pt'
vocab_file = f"{meta['resources_dir']}/source/vocab.json"
model_config_file = f"{meta['resources_dir']}/source/args.json"
input_file = mu.read(input)


def test_integration_embedding(run_component, tmp_path):
    output_embedding_file = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    run_component(
        [
            "--input",
            input,
            "--modality",
            "rna",
            "--model",
            model_file,
            "--model_vocab",
            vocab_file,
            "--model_config",
            model_config_file,
            "--dsbn",
            "True",
            "--obs_batch_label",
            "sample",
            "--obsm_gene_tokens",
            "gene_id_tokens",
            "--obsm_tokenized_values",
            "values_tokenized",
            "--obsm_padding_mask",
            "padding_mask",
            "--output",
            output_embedding_file,
        ]
    )

    # Read output file
    output_mdata = mu.read(output_embedding_file)
    output_adata = output_mdata.mod["rna"]

    # check that embedding obs is present
    assert (
        "X_scGPT" in output_adata.obsm.keys()
    ), "X_scGPT is not present in anndata obsm keys"

    # check embedding size
    assert (
        output_adata.obsm["X_scGPT"].shape[1] == 512
    ), "Embedding size does not equal 512"

    # check embedding value range
    assert not all(
        np.isnan(output_adata.obsm["X_scGPT"][0])
    ), "Embedding values are nan"
    assert all(
        [all(i > -1) & all(i < 1) for i in output_adata.obsm["X_scGPT"]]
    ), "Range of embedding values is outside of [-1, 1]"

    # Run embeddings without dsbn
    output_embedding_file_without_dsbn = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    run_component(
        [
            "--input",
            input,
            "--modality",
            "rna",
            "--model",
            model_file,
            "--model_vocab",
            vocab_file,
            "--model_config",
            model_config_file,
            "--dsbn",
            "False",
            "--obsm_gene_tokens",
            "gene_id_tokens",
            "--obsm_tokenized_values",
            "values_tokenized",
            "--obsm_padding_mask",
            "padding_mask",
            "--output",
            output_embedding_file_without_dsbn,
        ]
    )

    # Read output file
    output_mdata_no_dsbn = mu.read(output_embedding_file_without_dsbn)
    output_adata_no_dsbn = output_mdata_no_dsbn.mod["rna"]

    # Assert that embeddings without dsbn are different
    assert not (
        output_adata.obsm["X_scGPT"] == output_adata_no_dsbn.obsm["X_scGPT"]
    ).all(), "Embeddings with and without dsbn are the same"


def test_integration_embedding_dsbn_without_batch_labels(run_component, tmp_path):
    output_embedding_file = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    args = [
        "--input",
        input,
        "--modality",
        "rna",
        "--model",
        model_file,
        "--model_vocab",
        vocab_file,
        "--model_config",
        model_config_file,
        "--dsbn",
        "True",
        "--obsm_gene_tokens",
        "gene_id_tokens",
        "--obsm_tokenized_values",
        "values_tokenized",
        "--obsm_padding_mask",
        "padding_mask",
        "--output",
        output_embedding_file,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: When dsbn is set to True, you are required to provide batch labels \(input_obs_batch_labels\)\.",
        err.value.stdout.decode("utf-8"),
    )


def test_integration_embedding_non_existing_keys(run_component, tmp_path):
    output_embedding_file = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    # Test for non-existing gene names key
    args_1 = [
        "--input",
        input,
        "--modality",
        "rna",
        "--model",
        model_file,
        "--model_vocab",
        vocab_file,
        "--model_config",
        model_config_file,
        "--dsbn",
        "True",
        "--obs_batch_label",
        "sample",
        "--var_gene_names",
        "dummy_gene_name_key",
        "--obsm_gene_tokens",
        "gene_id_tokens",
        "--obsm_tokenized_values",
        "values_tokenized",
        "--obsm_padding_mask",
        "padding_mask",
        "--output",
        output_embedding_file,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_1)
    assert re.search(
        r"KeyError: \'dummy_gene_name_key\'", err.value.stdout.decode("utf-8")
    )

    # Test for non-existing batch label key
    args_2 = [
        "--input",
        input,
        "--modality",
        "rna",
        "--model",
        model_file,
        "--model_vocab",
        vocab_file,
        "--model_config",
        model_config_file,
        "--dsbn",
        "True",
        "--obs_batch_label",
        "dummy_batch_label_key",
        "--obsm_gene_tokens",
        "gene_id_tokens",
        "--obsm_tokenized_values",
        "values_tokenized",
        "--obsm_padding_mask",
        "padding_mask",
        "--output",
        output_embedding_file,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_2)
    assert re.search(
        r"KeyError: \'dummy_batch_label_key\'", err.value.stdout.decode("utf-8")
    )

    # Test for non-existing tokenized values key
    args_3 = [
        "--input",
        input,
        "--modality",
        "rna",
        "--model",
        model_file,
        "--model_vocab",
        vocab_file,
        "--model_config",
        model_config_file,
        "--dsbn",
        "True",
        "--obs_batch_label",
        "sample",
        "--obsm_gene_tokens",
        "gene_id_tokens",
        "--obsm_tokenized_values",
        "dummy_values_tokenized",
        "--obsm_padding_mask",
        "padding_mask",
        "--output",
        output_embedding_file,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_3)
    assert re.search(
        r'KeyError: "The parameter \'dummy_values_tokenized\' provided for \'--obsm_tokenized_values\' could not be found in adata.obsm"',
        err.value.stdout.decode("utf-8"),
    )


def test_finetuned_model(run_component, tmp_path):
    output_embedding_file = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    run_component(
        [
            "--input",
            input,
            "--modality",
            "rna",
            "--model",
            ft_model_file,
            "--model_vocab",
            vocab_file,
            "--model_config",
            model_config_file,
            "--dsbn",
            "True",
            "--obs_batch_label",
            "sample",
            "--obsm_gene_tokens",
            "gene_id_tokens",
            "--obsm_tokenized_values",
            "values_tokenized",
            "--obsm_padding_mask",
            "padding_mask",
            "--finetuned_checkpoints_key",
            "model_state_dict",
            "--output",
            output_embedding_file,
        ]
    )

    # Read output file
    output_mdata = mu.read(output_embedding_file)
    output_adata = output_mdata.mod["rna"]

    # check that embedding obs is present
    assert (
        "X_scGPT" in output_adata.obsm.keys()
    ), "X_scGPT is not present in anndata obsm keys"

    # check embedding size
    assert (
        output_adata.obsm["X_scGPT"].shape[1] == 512
    ), "Embedding size does not equal 512"

    # check embedding value range
    assert not all(
        np.isnan(output_adata.obsm["X_scGPT"][0])
    ), "Embedding values are nan"
    assert all(
        [all(i > -1) & all(i < 1) for i in output_adata.obsm["X_scGPT"]]
    ), "Range of embedding values is outside of [-1, 1]"


def test_finetuned_model_architecture(run_component, tmp_path):
    output_embedding_file = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    args = [
        "--input",
        input,
        "--modality",
        "rna",
        "--model",
        ft_model_file,
        "--model_vocab",
        vocab_file,
        "--model_config",
        model_config_file,
        "--dsbn",
        "True",
        "--obs_batch_label",
        "sample",
        "--obsm_gene_tokens",
        "gene_id_tokens",
        "--obsm_tokenized_values",
        "values_tokenized",
        "--obsm_padding_mask",
        "padding_mask",
        "--finetuned_checkpoints_key",
        "dummy_checkpoints_key",
        "--output",
        output_embedding_file,
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: The key \'dummy_checkpoints_key\' provided for \'--finetuned_checkpoints_key\' could not be found in the provided --model file. The finetuned model file for cell type annotation requires valid keys for the checkpoints and the label mapper.",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
