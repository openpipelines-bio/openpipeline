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
    "obsm_gene_tokens": 'gene_id_tokens',
    "obsm_tokenized_values": 'values_tokenized',
    "obsm_padding_mask": 'padding_mask',
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
    "dropout": 0.2,
    "DSBN": True,
    "n_input_bins": 51,
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

all_gene_ids = adata.obsm[par["obsm_gene_tokens"]]
all_values = adata.obsm[par["obsm_tokenized_values"]]
padding_mask = adata.obsm[par["obsm_padding_mask"]]

# Fetch batch ids for domain-specific batch normalization
if par["DSBN"] and not par["obs_batch_label"]:
    raise ValueError("When DSBN is set to True, you are required to provide batch labels (input_obs_batch_labels).")
elif par["DSBN"] and par["obs_batch_label"]:
    logger.info("Fetching batch id's for domain-specific batch normalization")
    batch_id_cats = adata.obs[par["obs_batch_label"]].astype("category")
    batch_id_labels = batch_id_cats.cat.codes.values
    batch_ids = batch_id_labels.tolist()
    batch_ids = np.array(batch_ids)
    num_batch_types = len(set(batch_ids))
elif not par["DSBN"] and par["obs_batch_label"]:
    logger.info("Batch labels provided but DSBN is set to False. Batch labels will be ignored and no DSBN will be performed.")

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
    dropout=par["dropout"],
    pad_token=pad_token,
    pad_value=pad_value,
    nlayers_cls=3,  # only applicable for decoder-based operations
    n_cls=1,  # only applicable for decoder-based operations
    do_mvc=False,  # only applicable for decoder-based operations
    ecs_threshold=0.8,  # only applicable for decoder-based operations
    do_dab=False,  # only applicable for decoder-based operations
    use_batch_labels=False, # only applicable for decoder-based operations
    num_batch_labels=num_batch_types if par["DSBN"] else None,
    domain_spec_batchnorm=par["DSBN"],
    input_emb_style="continuous",  # scGPT default
    explicit_zero_prob=False,  #TODO: Parametrize when GPU-based machine types are supported
    use_fast_transformer=False,  #TODO: Parametrize when GPU-based machine types are supported
    # fast_transformer_backend="flash",  #TODO: Parametrize when GPU-based machine types are supported
    pre_norm=False  #TODO: Parametrize when GPU-based machine types are supported
    )

load_pretrained(
    model,
    torch.load(model_file, map_location=device),
    verbose=False
    )

# Embed tokenized data
logger.info("Converting tokenized input data to embeddings")
model.to(device)
model.eval()

# cell_embeddings = model.encode_batch(
#     torch.from_numpy(all_gene_ids),
#     torch.from_numpy(all_values).float(),
#     src_key_padding_mask=torch.from_numpy(padding_mask),
#     batch_size=par["batch_size"],
#     batch_labels=torch.from_numpy(batch_ids).long() if par["DSBN"] else None,
#     output_to_cpu=True,
#     time_step=0,
#     return_np=True
# )

# cell_embeddings = cell_embeddings / np.linalg.norm(
#     cell_embeddings, axis=1, keepdims=True
# )

obs_shape = len(adata.obs["sample"])
cell_embeddings = np.zeros((obs_shape, 512))

# Write output
logger.info("Writing output data")
adata.obsm[par["obsm_embeddings"]] = cell_embeddings
mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])