import sys
import numpy as np
import mudata as mu
import json
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.utils.util import load_pretrained
import torch

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_subset_tokenized.h5mu",
    "obsm_gene_tokens": "gene_id_tokens",
    "obsm_tokenized_values": "values_tokenized",
    "obsm_padding_mask": "padding_mask",
    "model": "resources_test/scgpt/source/best_model.pt",
    "model_config": "resources_test/scgpt/source/args.json",
    "model_vocab": "resources_test/scgpt/source/vocab.json",
    "output": "Kim2020_Lung_embedded.h5ad",
    "var_gene_names": "gene_name",
    "obs_batch_label": "sample",
    "obsm_embeddings": "X_scGPT",
    "pad_token": "<pad>",
    "pad_value": -2,
    "batch_size": 64,
    "modality": "rna",
    "dsbn": True,
    "n_input_bins": 51,
}
meta = {
    "resources_dir": "src/utils",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info(f"Setting device to {'cuda' if torch.cuda.is_available() else 'cpu'}")
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

logger.info("Reading in data")

# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

for k, v in {
    "--obsm_gene_tokens": par["obsm_gene_tokens"],
    "--obsm_tokenized_values": par["obsm_tokenized_values"],
    "--obsm_padding_mask": par["obsm_padding_mask"],
}.items():
    if v not in adata.obsm.keys():
        raise KeyError(
            f"The parameter '{v}' provided for '{k}' could not be found in adata.obsm"
        )

all_gene_ids = adata.obsm[par["obsm_gene_tokens"]]
all_values = adata.obsm[par["obsm_tokenized_values"]]
padding_mask = adata.obsm[par["obsm_padding_mask"]]

# Fetch batch ids for domain-specific batch normalization
if par["dsbn"] and not par["obs_batch_label"]:
    raise ValueError(
        "When dsbn is set to True, you are required to provide batch labels (input_obs_batch_labels)."
    )
elif par["dsbn"] and par["obs_batch_label"]:
    logger.info("Fetching batch id's for domain-specific batch normalization")
    batch_id_cats = adata.obs[par["obs_batch_label"]].astype("category")
    batch_id_labels = batch_id_cats.cat.codes.values
    batch_ids = batch_id_labels.tolist()
    batch_ids = np.array(batch_ids)
    num_batch_types = len(set(batch_ids))
elif not par["dsbn"] and par["obs_batch_label"]:
    logger.info(
        "Batch labels provided but dsbn is set to False. Batch labels will be ignored and no dsbn will be performed."
    )

# Set padding specs
logger.info("Setting padding specs")
pad_token = par["pad_token"]
pad_value = par["pad_value"]
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Fetching gene names
logger.info("Fetching gene names")
if not par["var_gene_names"]:
    genes = adata.var.index.astype(str).tolist()
else:
    genes = adata.var[par["var_gene_names"]].astype(str).tolist()

# Model files
logger.info("Loading model, vocab and configs")
model_config_file = par["model_config"]
model_file = par["model"]
vocab_file = par["model_vocab"]

# Load vocab
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

vocab.set_default_index(vocab["<pad>"])
ntokens = len(vocab)
gene_ids = np.array(vocab(genes), dtype=int)

# Load model configs
with open(model_config_file, "r") as f:
    model_configs = json.load(f)
embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]

# Instantiate model
logger.info("Initializing transformer model")
model = TransformerModel(
    ntokens,
    d_model=embsize,
    nhead=nhead,
    d_hid=d_hid,
    nlayers=nlayers,
    vocab=vocab,
    dropout=0.5,  # scGPT default, only relevant for fine-tuning applications
    pad_token=pad_token,
    pad_value=pad_value,
    nlayers_cls=3,  # only applicable for decoder-based operations
    n_cls=1,  # only applicable for decoder-based operations
    do_mvc=False,  # only applicable for decoder-based operations
    ecs_threshold=0.8,  # only applicable for decoder-based operations
    do_dab=False,  # only applicable for decoder-based operations
    use_batch_labels=False,  # only applicable for decoder-based operations
    num_batch_labels=num_batch_types if par["dsbn"] else None,
    domain_spec_batchnorm=par["dsbn"],
    input_emb_style="continuous",  # scGPT default
    explicit_zero_prob=False,  # TODO: Parametrize when GPU-based machine types are supported
    use_fast_transformer=False,  # TODO: Parametrize when GPU-based machine types are supported
    # fast_transformer_backend="flash",  #TODO: Parametrize when GPU-based machine types are supported
    pre_norm=False,  # TODO: Parametrize when GPU-based machine types are supported
)


logger.info("Loading model")
model_file = par["model"]
model_dict = torch.load(model_file, map_location=device)

# Ensure the provided model has the correct architecture
finetuned_checkpoints_key = par.get("finetuned_checkpoints_key")
if finetuned_checkpoints_key:
    try:
        model_dict = model_dict[finetuned_checkpoints_key]
    except KeyError as e:
        raise ValueError(
            f"The key '{finetuned_checkpoints_key}' provided for '--finetuned_checkpoints_key' could not be found in the provided --model file. The finetuned model file for cell type annotation requires valid keys for the checkpoints and the label mapper."
        ) from e

# Load model
load_pretrained(model, model_dict, verbose=False)

# Embed tokenized data
logger.info("Converting tokenized input data to embeddings")
model.to(device)
model.eval()

cell_embeddings = model.encode_batch(
    torch.from_numpy(all_gene_ids),
    torch.from_numpy(all_values).float(),
    src_key_padding_mask=torch.from_numpy(padding_mask),
    batch_size=par["batch_size"],
    batch_labels=torch.from_numpy(batch_ids).long() if par["dsbn"] else None,
    output_to_cpu=True,
    time_step=0,
    return_np=True,
)

cell_embeddings = cell_embeddings / np.linalg.norm(
    cell_embeddings, axis=1, keepdims=True
)

# Write output
logger.info("Writing output data")
adata.obsm[par["obsm_embeddings"]] = cell_embeddings
mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])
