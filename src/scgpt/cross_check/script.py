import anndata as ad
import numpy as np
from pathlib import Path
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)

## VIASH START
par = {
    "input": "src/scgpt/test_resources/Kim2020_Lung.h5ad",
    "output": "src/scgpt/test_resources/Kim2020_Lung_crosschecked.h5ad",
    "ori_batch_layer_name": "sample",
    "batch_id_layer": "batch_id",
    "gene_name_layer": "gene_name",
    "pad_token": "<pad>",
    "load_model_vocab": True,
    "model_dir": "src/scgpt/model"
}
## VIASH END

# Read in data
input_adata = ad.read_h5ad(par["input"])
adata = input_adata.copy()

pad_token = par["pad_token"]
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Make batch a category column
adata.obs["str_batch"] = adata.obs[par["ori_batch_layer_name"]].astype(str)
batch_id_labels = adata.obs["str_batch"].astype("category").cat.codes.values
adata.obs[par["batch_id_layer"]] = batch_id_labels
adata.var[par["gene_name_layer"]] = adata.var.index.tolist()

# Cross-check genes with pre-trained model
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

adata.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in adata.var[par["gene_name_layer"]]
    ]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]

adata.write_h5ad(par["output"])
