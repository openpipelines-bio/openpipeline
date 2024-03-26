import sys
import os
import mudata as mu
import numpy as np
from scgpt.tokenizer.gene_tokenizer import GeneVocab

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_subset.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "ori_batch_layer_name": "sample",
    "batch_id_layer": "batch_id",
    "gene_name_layer": "gene_name",
    "pad_token": "<pad>",
    "vocab_file": "resources_test/scgpt/source/vocab.json"
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

# Read in data
logger.info(f"Reading {par['input']}")
mudata = mu.read_h5mu(par["input"])
adata = mudata.mod[par["modality"]].copy()

pad_token = par["pad_token"]
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Make batch a category column
logger.info(f"Converting batch column '{par['ori_batch_layer_name']}' to category")

if par["ori_batch_layer_name"] not in adata.obs.columns:
    raise ValueError(f"Batch column '{par['ori_batch_layer_name']}' not found in .mod['{par['modality']}'].obs.")

adata.obs["str_batch"] = adata.obs[par["ori_batch_layer_name"]].astype(str).tolist()
batch_id_labels = adata.obs["str_batch"].astype("category").cat.codes.values
adata.obs[par["batch_id_layer"]] = batch_id_labels.astype(str).tolist()
adata.var[par["gene_name_layer"]] = adata.var.index.astype(str).tolist()

# Cross-check genes with pre-trained model
logger.info("Cross-checking genes with pre-trained model")
genes = adata.var[par["gene_name_layer"]].tolist()

logger.info(f"Loading model vocab from {par['vocab_file']}")

vocab_file = par["vocab_file"]
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

logger.info("Filtering genes based on model vocab")
adata.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in adata.var[par["gene_name_layer"]]
    ]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]

mudata.mod[par["modality"]] = adata

logger.info(f"Writing to {par['output']}")
mudata.write_h5mu(par["output"])
