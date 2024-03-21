import mudata as mu
import numpy as np
from pathlib import Path
import torch
from scipy.sparse import issparse
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab


import pytest
import sys
import mudata as mu
import torch
import numpy as np
from pathlib import Path
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor

## VIASH START
meta = {
    "resources_dir": "resources_test",
    "executable": "./target/docker/scgpt/integration_embedding/integration_embedding",
    "temp_dir": "tmp",
    "config": "./target/docker/scgpt/integration_embedding/.config.vsh.yaml"
}
## VIASH END

input = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung.h5mu"
model_dir = f"{meta['resources_dir']}/scgpt/source/"
input_file = mu.read(input)

## START TEMPORARY WORKAROUND (until all scgpt modules are implemented)
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
model_dir = Path(model_dir)
vocab_file = model_dir / "vocab.json"
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

#input_file.mod["rna"] = adata
#########

all_counts = (
    adata.layers["X_binned"].A
    if issparse(adata.layers["X_binned"])
    else adata.layers["X_binned"]
)

# Fetch gene names and look up tokens in vocab
vocab.set_default_index(vocab["<pad>"])
ntokens = len(vocab)
genes = adata.var["gene_name"].tolist()
gene_ids = np.array(vocab(genes), dtype=int)

# Fetch number of subset hvg
n_hvg = adata.var.shape[0]

# Tokenize and pad data
tokenized_data = tokenize_and_pad_batch(
    all_counts,
    gene_ids,
    max_len=n_hvg+1,
    vocab=vocab,
    pad_token=pad_token,
    pad_value=-2,
    append_cls=True,  # append <cls> token at the beginning,
    include_zero_gene=False,
    return_pt=True,
    mod_type=None,
    vocab_mod=None
    )

all_gene_ids, all_values = tokenized_data["genes"], tokenized_data["values"]
padding_mask = all_gene_ids.eq(vocab[pad_token])

## END OF TEMPORARY WORKAROUND

def test_integration_embedding(run_component, tmp_path):

    input_gene_ids = tmp_path / "gene_ids.pt"
    input_values = tmp_path / "values.pt"
    input_padding_mask = tmp_path / "padding_mask.pt"

    torch.save(all_gene_ids, input_gene_ids)
    torch.save(all_values, input_values)
    torch.save(padding_mask, input_padding_mask)

    input_preprocessed = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_preprocessed.h5mu"
    input_file.write(input_preprocessed)

    output_embedding_file = tmp_path / "Kim2020_Lung_embedded.h5mu"

    run_component([
        "--input", input_preprocessed,
        "--modality", "rna",
        "--input_gene_ids", input_gene_ids,
        "--input_values", input_values,
        "--input_padding_mask", input_padding_mask,
        "--pad_token", "<pad>",
        "--pad_value", "-2",
        "--dropout", 0.2,
        "--DSBN", True,
        "--GEPC", True,
        "--n_input_bins", 51,
        "--ecs_threshold", 0.8,
        "--explicit_zero_prob", True,
        "--use_fast_transformer", False,
        "--pre_norm", False,
        "--device", "cpu",
        "--batch_size", 64,
        "--input_layer", "X_binned",
        "--gene_name_layer", "gene_name",
        "--batch_id_layer", "batch_id",
        "--output", output_embedding_file,
        "--model_dir", model_dir
    ])

    # check that embedding obs is present
    assert 'X_scGPT' in input_preprocessed.obsm.keys(), "X_scGPT is not present in anndata obsm keys"

    # check dimensions
    assert input_preprocessed.obsm["X_scGPT"].shape[1] == 512, "Embedding size does not equal 512"
    assert input_preprocessed.obsm["X_scGPT"].shape[0] == all_gene_ids.shape[0], "Embedding dimensions don't match input adata dimension"

    # check values
    assert not all(np.isnan(adata.obsm["X_scGPT"][0])), "Embedding values are nan"


# if __name__ == '__main__':
#     sys.exit(pytest.main([__file__]))
