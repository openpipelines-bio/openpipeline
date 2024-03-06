import anndata as ad
import numpy as np
from pathlib import Path
import json
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.utils.util import load_pretrained
import torch
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)

## VIASH START
par = {
    "input": "src/scgpt/data/Kim2020_Lung_preprocessed.h5ad",
    "input_gene_ids": 'src/scgpt/data/Kim2020_Lung_gene_ids.pt',
    "input_values": 'src/scgpt/data/Kim2020_Lung_values.pt',
    "input_padding_mask": 'src/scgpt/data/Kim2020_Lung_padding_mask.pt',
    "output": "src/scgpt/data/Kim2020_Lung_embedded.h5ad",
    "gene_name_layer": "gene_name",
    "batch_id_layer": "batch_id",
    "embedding_layer": "X_scGPT",
    "pad_token": "<pad>",
    "pad_value": -2,
    "load_model_vocab": False,
    "load_model_configs": False,
    "model_dir": "src/scgpt/model",
    "layer_size": 128,
    "nhead": 4,
    "nlayers": 4,
    "dropout": 0.2,
    "GEPC": True,
    "DSBN": True,
    "n_input_bins": 51,
    "ecs_threshold": 0.8,
    "explicit_zero_prob": True,
    "use_fast_transformer": False,
    "pre_norm": False,
    "device": "cpu",  # torch.device type
    "batch_size": 64
}
## VIASH END

# Read in data
input_adata = ad.read_h5ad(par["input"])
adata = input_adata.copy()
all_gene_ids = torch.load(par["input_gene_ids"])
all_values = torch.load(par["input_values"])
padding_mask = torch.load(par["input_padding_mask"])

# Fetch batch ids
batch_ids = adata.obs[par["batch_id_layer"]].tolist()
batch_ids = np.array(batch_ids)
num_batch_types = len(set(batch_ids))

# Set padding specs
pad_token = par["pad_token"]
pad_value = par["pad_value"]
special_tokens = [pad_token, "<cls>", "<eoc>"]
genes = adata.var[par["gene_name_layer"]].tolist()

# Load vocab
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

# Set model configs
model_dir = Path(par["model_dir"])
model_config_file = model_dir / "args.json"
model_file = model_dir / "best_model.pt"

if par["load_model_configs"]:
    with open(model_config_file, "r") as f:
        model_configs = json.load(f)
    embsize = model_configs["embsize"]
    nhead = model_configs["nheads"]
    d_hid = model_configs["d_hid"]
    nlayers = model_configs["nlayers"]
    n_layers_cls = model_configs["n_layers_cls"]
else:
    embsize = par["layer_size"]
    nhead = par["nhead"]
    nlayers = par["nlayers"]
    d_hid = par["layer_size"]

# Load model
model = TransformerModel(
    ntokens,
    embsize,
    nhead,
    d_hid,
    nlayers,
    vocab=vocab,
    dropout=par["dropout"],
    pad_token=pad_token,
    pad_value=pad_value,
    do_mvc=par["GEPC"],
    do_dab=True,
    use_batch_labels=True,
    num_batch_labels=num_batch_types,
    domain_spec_batchnorm=par["DSBN"],
    n_input_bins=par["n_input_bins"],
    ecs_threshold=par["ecs_threshold"],
    explicit_zero_prob=par["explicit_zero_prob"],
    use_fast_transformer=par["use_fast_transformer"],
    pre_norm=par["pre_norm"],
)


if par["device"] == "cuda" and not torch.cuda.is_avilable():
    print("WARNING: CUDA is not available. Using CPU instead.")    
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device = torch.device(par["device"])

load_pretrained(model, torch.load(model_file, map_location=device), verbose=False)

model.to(device)
model.eval()

# Embed tokenized data
cell_embeddings = model.encode_batch(
    all_gene_ids,
    all_values.float(),
    src_key_padding_mask=padding_mask,
    batch_size=par["batch_size"],
    batch_labels=torch.from_numpy(batch_ids).long() if par["DSBN"] else None,
    time_step=0,
    return_np=True
)
cell_embeddings = cell_embeddings / np.linalg.norm(
    cell_embeddings, axis=1, keepdims=True
)

adata.obsm[par["embedding_layer"]] = cell_embeddings
adata.write_h5ad(par["output"])
