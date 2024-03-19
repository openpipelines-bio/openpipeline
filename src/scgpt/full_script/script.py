import anndata as ad
import numpy as np
from pathlib import Path
import json
from scipy.sparse import issparse
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.preprocess import Preprocessor
from scgpt.utils.util import load_pretrained
import torch
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)

## VIASH START
par = {
    "input": "src/scgpt/data/filtered_ms_adata.h5ad",
    "output": "src/scgpt/data/filtered_ms_adata_embedded.h5ad",
    "modality": "rna",
    "input_layer": "X_binned",
    "gene_name_layer": "gene_name",
    "batch_id_layer": "batch_id",
    "pad_token": "<pad>",
    "pad_value": -2,
    "n_hvg": 1200,
    "DSBN": True,
    "load_model": None,
    "model_dir": "/Users/dorienroosen/code/openpipeline/src/scgpt/model",
    "device": "cpu",  # torch.device type
    "nhead": 4,
    "nlayers": 4,
    "dropout": 0.2,
    "GEPC": True,
    "n_input_bins": 51,
    "ecs_threshold": 0.8,
    "explicit_zero_prob": True,
    "use_fast_transformer": False,
    "pre_norm": False,
    "filter_gene_by_counts": 3,
    "filter_cell_by_counts": False,
    "normalize_total": 1e4,
    "normalized_total_layer": "X_normed",
    "log1p_layer": "X_log1p",
    "subset_hvg": False,
    "data_is_raw": False,
    "per_seq_batch_sample": False,
    "layer_size": 128,
    "batch_size": 64
}
## VIASH END

# Read in data
input_adata = ad.read_h5ad(par["input"])
adata = input_adata.copy()

# Read in model files
model_dir = Path(par["model_dir"])
model_config_file = model_dir / "args.json"
vocab_file = model_dir / "vocab.json"
model_file = model_dir / "best_model.pt"

# Define various layer names
input_layer_key = par["input_layer"]
gene_name_layer = par["gene_name_layer"]
batch_id_layer = par["batch_id_layer"]
# Set model configs
if par["load_model"] is not None:
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

# Set padding specs
pad_token = par["pad_token"]
special_tokens = [pad_token, "<cls>", "<eoc>"]
pad_value = -2

# Preprocess data
preprocessor = Preprocessor(
    use_key="X",  # the key in adata.layers to use as raw data
    filter_gene_by_counts=par["filter_gene_by_counts"],  # step 1
    filter_cell_by_counts=par["filter_cell_by_counts"],  # step 2
    normalize_total=par["normalize_total"],  # 3. whether to normalize the raw data and to what sum
    result_normed_key=par["normalized_total_layer"],  # the key in adata.layers to store the normalized data
    log1p=par["data_is_raw"],  # 4. whether to log1p the normalized data
    result_log1p_key=par["log1p_layer"],
    subset_hvg=par["n_hvg"],  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor="seurat_v3" if par["data_is_raw"] else "cell_ranger",
    binning=par["n_input_bins"],  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
    )

preprocessor(adata, batch_key="str_batch")
if par["per_seq_batch_sample"]:
    # sort the adata by batch_id in advance
    adata_sorted = adata[adata.obs["batch_id"].argsort()].copy()

# Fetch counts
all_counts = (
    adata.layers[input_layer_key].A
    if issparse(adata.layers[input_layer_key])
    else adata.layers[input_layer_key]
)

# Fetch gene names and look up tokens in vocab
genes = adata.var[gene_name_layer].tolist()

if par["load_model"] is None:
    vocab = Vocab(
        VocabPybind(genes + special_tokens, None)
    )  # bidirectional lookup [gene <-> int]
elif par["load_model"] is not None:
    vocab = GeneVocab.from_file(vocab_file)
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)

vocab.set_default_index(vocab["<pad>"])
ntokens = len(vocab)
gene_ids = np.array(vocab(genes), dtype=int)

# Fetch batch ids
batch_ids = adata.obs[batch_id_layer].tolist()
batch_ids = np.array(batch_ids)
num_batch_types = len(set(batch_ids))

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

# Load model
model = TransformerModel(
    ntokens,
    embsize,
    par["nhead"],
    d_hid,
    par["nlayers"],
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

adata.obsm["X_scGPT"] = cell_embeddings
adata.write_h5ad(par["output"])
