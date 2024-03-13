import anndata as ad
import numpy as np
from pathlib import Path
import torch
from scipy.sparse import issparse
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)

## VIASH START
par = {
    "input": "src/scgpt/test_resources/Kim2020_Lung_preprocessed.h5ad",
    "output_gene_ids": 'src/scgpt/test_resources/Kim2020_Lung_gene_ids.pt',
    "output_values": 'src/scgpt/test_resources/Kim2020_Lung_values.pt',
    "output_padding_mask": 'src/scgpt/test_resources/Kim2020_Lung_padding_mask.pt',
    "pad_token": "<pad>",
    "pad_value": -2,
    "input_layer": "X_binned",
    "gene_name_layer": "gene_name",
    "load_model_vocab": False,
    "model_dir": "src/scgpt/model",
    "n_hvg": 1200
}
## VIASH END

# Read in data
input_adata = ad.read_h5ad(par["input"])
adata = input_adata.copy()

# Set padding specs
pad_token = par["pad_token"]
special_tokens = [pad_token, "<cls>", "<eoc>"]
pad_value = -2

# Fetch counts
all_counts = (
    adata.layers[par["input_layer"]].A
    if issparse(adata.layers[par["input_layer"]])
    else adata.layers[par["input_layer"]]
)

# Fetch gene names and look up tokens in vocab
genes = adata.var[par["gene_name_layer"]].tolist()

if par["load_model_vocab"]:
    model_dir = Path(par["model_dir"])
    vocab_file = model_dir / "vocab.json"
    vocab = GeneVocab.from_file(vocab_file)
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
else:
    # bidirectional lookup [gene <-> int]
    vocab = Vocab(
        VocabPybind(genes + special_tokens, None)
    )

vocab.set_default_index(vocab["<pad>"])
ntokens = len(vocab)
gene_ids = np.array(vocab(genes), dtype=int)

# Tokenize and pad data
tokenized_data = tokenize_and_pad_batch(
    all_counts,
    gene_ids,
    max_len=par["n_hvg"]+1,
    vocab=vocab,
    pad_token=pad_token,
    pad_value=pad_value,
    append_cls=True,  # append <cls> token at the beginning,
    include_zero_gene=False,
    return_pt=True,
    mod_type=None,
    vocab_mod=None
    )

all_gene_ids, all_values = tokenized_data["genes"], tokenized_data["values"]
padding_mask = all_gene_ids.eq(vocab[pad_token])

torch.save(all_gene_ids, par["output_gene_ids"])
torch.save(all_values, par["output_values"])
torch.save(padding_mask, par["output_padding_mask"])