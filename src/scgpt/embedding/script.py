import numpy as np
import mudata as mu
import json
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.utils.util import load_pretrained
import torch

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_tokenized.h5mu",
    "input_obsm_gene_tokens": 'gene_id_tokens',
    "input_obsm_tokenized_values": 'values_tokenized',
    "input_obsm_padding_mask": 'padding_mask',
    "model": "resources_test/scgpt/source/best_model.pt",
    "model_config": "resources_test/scgpt/source/args.json",
    "model_vocab": "resources_test/scgpt/source/vocab.json",
    "output": "Kim2020_Lung_embedded.h5ad",
    "input_var_gene_names": "gene_name",
    "input_obs_batch_label": "sample",
    "embedding_layer_key": "X_scGPT",
    "pad_token": "<pad>",
    "pad_value": -2,
    "modality": "rna",
    "dropout": 0.2,
    "GEPC": True,
    "DSBN": True,
    "n_input_bins": 51,
    "ecs_threshold": 0.8,
    "explicit_zero_prob": True,
    "use_fast_transformer": False,
    "pre_norm": False,
    "batch_size": 64,
    "output_compression": None
}
## VIASH END

# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()


#TODO: Optionally set device to use gpu, to be implemented when gpu-based machine types are supported
# if par["device"] == "cuda" and not torch.cuda.is_available():
#     print("WARNING: CUDA is not available. Using CPU instead.")
#     device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# device = torch.device(par["device"])

logger.info("Setting device to use cpu")

device = torch.device("cpu")

logger.info("Reading in data")

# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

all_gene_ids = adata.obsm[par["input_obsm_gene_tokens"]]
all_values = adata.obsm[par["input_obsm_tokenized_values"]]
padding_mask = adata.obsm[par["input_obsm_padding_mask"]]

# Fetch batch ids for domain-specific batch normalization
if par["DSBN"]:
    if not par["input_obs_batch_label"]:
        raise ValueError("When DSBN is set to True, you are required to provide batch labels (input_obs_batch_labels).")
    else:
        batch_id_cats = adata.obs[par["input_obs_batch_label"]].astype("category")
        batch_id_labels = batch_id_cats.cat.codes.values
        batch_ids = batch_id_labels.tolist()
        batch_ids = np.array(batch_ids)
        num_batch_types = len(set(batch_ids))

# Set padding specs
pad_token = par["pad_token"]
pad_value = par["pad_value"]
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Fetching gene names
if not par["input_var_gene_names"]:
    genes = adata.var.index.astype(str).tolist()
else:
    genes = adata.var[par["input_var_gene_names"]].astype(str).tolist()

logger.info("Loading model, vocab and configs")
# Model files
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

logger.info("Initializing transformer model")

# Instantiate model
model = TransformerModel(
    ntokens,
    d_model=embsize,
    nhead=nhead,
    d_hid=d_hid,
    nlayers=nlayers,
    vocab=vocab,
    dropout=par["dropout"],
    pad_token=pad_token,
    pad_value=pad_value,
    nlayers_cls=3,  #TODO: parametrize for decoder-based operations
    n_cls=1,  #TODO: parametrize for decoder-based operations
    do_mvc=True,  #TODO: parametrize for decoder-based operations
    ecs_threshold=0.8,  #TODO: parametrize for decoder-based operations
    do_dab=True, #TODO: parametrize for decoder-based operations
    use_batch_labels=True,
    num_batch_labels=num_batch_types if par["DSBN"] else None,
    domain_spec_batchnorm=par["DSBN"],
    input_emb_style="continuous",  # scGPT default
    # n_input_bins=par["n_input_bins"],  # only applies when input_emb_style is "category"
    explicit_zero_prob=True,  #TODO: Parametrize when GPU-based machine types are supported
    use_fast_transformer=False,  #TODO: Parametrize when GPU-based machine types are supported
    fast_transformer_backend="flash",  #TODO: Parametrize when GPU-based machine types are supported
    pre_norm=False  #TODO: Parametrize when GPU-based machine types are supported
    )

load_pretrained(
    model,
    torch.load(model_file, map_location=device),
    verbose=False
    )

model.to(device)
model.eval()

logger.info("Converting tokenized input data to embeddings")
# Embed tokenized data
cell_embeddings = model.encode_batch(
    torch.from_numpy(all_gene_ids),
    torch.from_numpy(all_values).float(),
    src_key_padding_mask=torch.from_numpy(padding_mask),
    batch_size=par["batch_size"],
    batch_labels=torch.from_numpy(batch_ids).long() if par["DSBN"] else None,
    output_to_cpu=True,
    time_step=0,
    return_np=True
)

cell_embeddings = cell_embeddings / np.linalg.norm(
    cell_embeddings, axis=1, keepdims=True
)

logger.info("Writing output data")
# Write output
adata.obsm[par["embedding_layer_key"]] = cell_embeddings
mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])