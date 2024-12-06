import pytest
import sys
import mudata as mu
from scgpt.tokenizer.gene_tokenizer import GeneVocab

## VIASH START
meta = {
    "resources_dir": "resources_test/scgpt",
    "executable": "./target/docker/scgpt/integration_pad_tokenize/integration_pad_tokenize",
    "temp_dir": "tmp",
    "config": "./target/docker/scgpt/integration_pad_tokenize/.config.vsh.yaml",
}
## VIASH END

input_file = (
    f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_subset_binned.h5mu"
)
vocab_file = f"{meta['resources_dir']}/scgpt/source/vocab.json"
vocab = GeneVocab.from_file(vocab_file)


def test_integration_pad_tokenize(run_component, tmp_path):
    output = tmp_path / "Kim2020_Lung_tokenized.h5mu"

    run_component([
        "--input", input_file,
        "--output", output,
        "--modality", "rna",
        "--var_input", "scgpt_cross_checked_genes",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--obsm_padding_mask", "padding_mask",
        "--pad_token", "<pad>",
        "--pad_value", "-2",
        "--input_obsm_binned_counts", "binned_counts",
        "--model_vocab", vocab_file
    ])

    output_file = mu.read(output)
    output_adata = output_file.mod["rna"]

    gene_ids = output_adata.obsm["gene_id_tokens"]
    values = output_adata.obsm["values_tokenized"]
    padding_mask = output_adata.obsm["padding_mask"]

    # check output dimensions
    ## nr of genes that are tokenized
    assert (
        gene_ids.shape[1] <= output_adata.var.shape[0] + 1
    ), "gene_ids shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"
    assert (
        values.shape[1] <= output_adata.var.shape[0] + 1
    ), "values shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"
    assert (
        padding_mask.shape[1] <= output_adata.var.shape[0] + 1
    ), "padding_mask shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"

    ## equal size of output tensors
    assert (
        gene_ids.shape == values.shape
    ), "gene_ids shape[1] does not match values shape[1]"
    assert (
        gene_ids.shape == padding_mask.shape
    ), "gene_ids shape[1] does not match padding_mask shape[1]"

    ## check values of output tensors
    assert gene_ids.dtype == "int64", "tokenized gene_ids are not integers"
    assert (gene_ids > 0).all(), "not all gene id tokens are higher than 0"

    assert values.dtype == "float32", "tokenized values are not floats"
    assert (values >= -2).all(), "not all tokenized values are higher than/equal to -2"

    assert padding_mask.dtype == bool, "padding mask is not boolean"

    ## check cls token
    assert (
        gene_ids[:, 0] == vocab["<cls>"]
    ).all(), (
        "cls token was not correctly appended at the beginning of the gene_ids tensor"
    )
    assert (
        values[:, 0] == 0
    ).all(), (
        "cls token was not correctly appended at the beginning of the values tensors"
    )

    # check padding values
    masked_gene_ids = gene_ids[padding_mask]
    unmasked_gene_ids = gene_ids[~padding_mask]
    assert all(
        masked_gene_ids == vocab["<pad>"]
    ), "masked gene_ids contain non-pad tokens"
    assert all(
        unmasked_gene_ids != vocab["<pad>"]
    ), "unmasked gene_ids contain pad tokens"

    masked_values = values[padding_mask]
    unmasked_values = values[~padding_mask]
    assert all(masked_values == -2), "masked values contain non-pad values"
    assert all(unmasked_values != -2), "unmasked values contain pad values"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
