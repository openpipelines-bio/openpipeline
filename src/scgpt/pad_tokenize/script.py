import mudata as mu
import numpy as np
from scipy.sparse import issparse
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab


## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_preprocessed.h5mu",
    "model_vocab": "resources_test/scgpt/source/vocab.json",
    "output": "resources_test/scgpt/test_resources/Kim2020_Lung_tokenized.h5mu",
    "pad_token": "<pad>",
    "pad_value": -2,
    "modality": "rna",
    "input_layer": "X_binned",
    "max_seq_len": None,
    "var_gene_names": None,
    "obsm_gene_tokens": "gene_id_tokens",
    "obsm_tokenized_values": "values_tokenized",
    "obsm_padding_mask": "padding_mask",
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

logger.info("Reading in data")

# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

# Set padding specs
pad_token = par["pad_token"]
special_tokens = [pad_token, "<cls>", "<eoc>"]
pad_value = -2

logger.info("Fetching counts and gene names")
# Fetch counts
all_counts = (
    adata.layers[par["input_layer"]].A
    if issparse(adata.layers[par["input_layer"]])
    else adata.layers[par["input_layer"]]
)

# Fetching gene names
if not par["var_gene_names"]:
    genes = adata.var.index.astype(str).tolist()
else:
    genes = adata.var[par["var_gene_names"]].astype(str).tolist()

# Fetch gene names and look up tokens in vocab
logger.info("Reading in vocab and fetching gene tokens")
vocab_file = par["model_vocab"]
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

vocab.set_default_index(vocab["<pad>"])
ntokens = len(vocab)
gene_ids = np.array(vocab(genes), dtype=int)

# Fetch max seq len
if not par["max_seq_len"]:
    max_seq_len = adata.var.shape[0] + 1
else:
    max_seq_len = par["max_seq_len"]

# Tokenize and pad data
logger.info("Padding and tokenizing data")
tokenized_data = tokenize_and_pad_batch(
    all_counts,
    gene_ids,
    max_len=max_seq_len,
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

logger.info("Writing output data")
adata.obsm[par["obsm_gene_tokens"]] = all_gene_ids.numpy()
adata.obsm[par["obsm_tokenized_values"]] = all_values.numpy()
adata.obsm[par["obsm_padding_mask"]] = padding_mask.numpy()

mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])
