import pytest
import sys
import mudata as mu
import torch
import numpy as np
from pathlib import Path
from scipy.sparse import issparse
from scgpt.tokenizer import tokenize_and_pad_batch
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

input = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_subset.h5mu"
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

input_gene_id_path = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_gene_ids.pt"
input_values_path = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_values.pt"
input_padding_mask_path = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_padding_mask.pt"

torch.save(all_gene_ids, input_gene_id_path)
torch.save(all_values, input_values_path)
torch.save(padding_mask, input_padding_mask_path)

input_preprocessed = mu.MuData({'rna': adata})
input_preprocessed_path = f"{meta['resources_dir']}/Kim2020_Lung_preprocessed.h5mu"
input_preprocessed.write_h5mu(input_preprocessed_path)

## END OF TEMPORARY WORKAROUND


def test_integration_embedding(run_component, tmp_path):

    output_embedding_file = tmp_path / "Kim2020_Lung_embedded.h5mu"

    run_component([
        "--input", input_preprocessed_path,
        "--modality", "rna",
        "--model_dir", model_dir,
        "--input_gene_ids", input_gene_id_path,
        "--input_values", input_values_path,
        "--input_padding_mask", input_padding_mask_path,
        "--output", output_embedding_file,
    ])

    # check that embedding obs is present
    assert 'X_scGPT' in input_preprocessed.obsm.keys(), "X_scGPT is not present in anndata obsm keys"

    # check dimensions
    assert input_preprocessed.obsm["X_scGPT"].shape[1] == 512, "Embedding size does not equal 512"
    assert input_preprocessed.obsm["X_scGPT"].shape[0] == all_gene_ids.shape[0], "Embedding dimensions don't match input adata dimension"

    # check values are not nan
    assert not all(np.isnan(adata.obsm["X_scGPT"][0])), "Embedding values are nan"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
