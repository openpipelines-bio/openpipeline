import pytest
import sys
import mudata as mu
import numpy as np
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor

## VIASH START
meta = {
    "resources_dir": "resources_test/scgpt",
    "executable": "./target/docker/scgpt/integration_pad_tokenize/integration_pad_tokenize",
    "temp_dir": "tmp",
    "config": "./target/docker/scgpt/integration_pad_tokenize/.config.vsh.yaml"
}
## VIASH END

input = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_subset.h5mu"
vocab_file = f"{meta['resources_dir']}/scgpt/source/vocab.json"
input_file = mu.read(input)

## START TEMPORARY WORKAROUND DATA PREPROCESSING
#TODO: Remove this workaround once full scGPT preprocessing workflow is implemented
# Read in data
adata = input_file.mod["rna"]

# Set tokens for integration
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Make batch a category column
adata.obs["str_batch"] = adata.obs["sample"].astype(str)
batch_id_labels = adata.obs["str_batch"].astype("category").cat.codes.values
adata.obs["batch_id"] = batch_id_labels
adata.var["gene_name"] = adata.var.index.tolist()

# Load model vocab
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

# Cross-check genes with pre-trained model
genes = adata.var["gene_name"].tolist()
adata.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in adata.var["gene_name"]
    ]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]

# Preprocess data
preprocessor = Preprocessor(
    use_key="X",
    filter_gene_by_counts=3,
    filter_cell_by_counts=False,
    normalize_total=10000,
    result_normed_key="X_normed",
    log1p=True,
    result_log1p_key="X_log1p",
    subset_hvg=1200,
    hvg_flavor="seurat_v3",
    binning=51,
    result_binned_key="X_binned",
    )

preprocessor(adata, batch_key="str_batch")

# copy results to mudata
input_file.mod["rna"] = adata

## END TEMPORARY WORKAROUND DATA PREPROCESSING


def test_integration_pad_tokenize(run_component, tmp_path):
    output = tmp_path / "Kim2020_Lung_tokenized.h5mu"
    input_preprocessed = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_preprocessed.h5mu"
    input_file.write(input_preprocessed)

    run_component([
        "--input", input_preprocessed,
        "--output", output,
        "--modality", "rna",
        "--output_obsm_gene_tokens", "gene_id_tokens",
        "--output_obsm_tokenized_values", "values_tokenized",
        "--output_obsm_padding_mask", "padding_mask",
        "--pad_token", "<pad>",
        "--pad_value", "-2",
        "--input_layer", "X_binned",
        "--gene_name_layer", "gene_name",
        "--model_vocab", vocab_file
    ])

    output_file = mu.read(output)
    output_adata = output_file.mod["rna"]

    gene_ids = output_adata.obsm["output_obsm_gene_tokens"]
    values = output_adata.obsm["output_obsm_tokenized_values"]
    padding_mask = output_adata.obsm["output_obsm_padding_mask"]

    # check output dimensions
    ## nr of genes that are tokenized 
    assert gene_ids.shape[1] <= output_adata.var.shape[0] + 1, "gene_ids shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"
    assert values.shape[1] <= output_adata.var.shape[0] + 1, "values shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"
    assert padding_mask.shape[1] <= output_adata.var.shape[0] + 1, "padding_mask shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"

    ## equal size of output tensors
    assert gene_ids.shape == values.shape, "gene_ids shape[1] does not match values shape[1]"
    assert gene_ids.shape == padding_mask.shape, "gene_ids shape[1] does not match padding_mask shape[1]"

    ## check values of output tensors
    assert gene_ids.dtype == "int64", "tokenized gene_ids are not integers"
    assert (gene_ids > 0).all(), "not all gene id tokens are higher than 0"

    assert values.dtype == "float32", "tokenized values are not floats"
    assert (values >= -2).all(), "not all tokenized values are higher than/equal to -2"

    assert padding_mask.dtype == bool, "padding mask is not boolean"

    ## check cls token
    assert (gene_ids[:, 0] == vocab["<cls>"]).all(), "cls token was not correctly appended at the beginning of the gene_ids tensor"
    assert (values[:, 0] == 0).all(), "cls token was not correctly appended at the beginning of the values tensors"

    # check padding values
    masked_gene_ids = gene_ids[padding_mask]
    unmasked_gene_ids = gene_ids[~padding_mask]
    assert all(masked_gene_ids == vocab["<pad>"]), "masked gene_ids contain non-pad tokens"
    assert all(unmasked_gene_ids != vocab["<pad>"]), "unmasked gene_ids contain pad tokens"

    masked_values = values[padding_mask]
    unmasked_values = values[~padding_mask]
    assert all(masked_values == -2), "masked values contain non-pad values"
    assert all(unmasked_values != -2), "unmasked values contain pad values"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
