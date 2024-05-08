import mudata as mu
import numpy as np
from scgpt.tokenizer.gene_tokenizer import GeneVocab

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_subset.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "input_var_gene_names": None,
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

# Fetching gene names
if not par["var_gene_names"]:
    genes = adata.var.index.astype(str).tolist()
elif par["var_gene_names"] not in adata.var.columns:
    raise ValueError(f"Gene name column '{par['var_gene_names']}' not found in .mod['{par['modality']}'].obs.")
else: 
    genes = adata.var[par["var_gene_names"]].astype(str).tolist()

# Cross-check genes with pre-trained model
logger.info(f"Loading model vocab from {par['vocab_file']}")
vocab_file = par["vocab_file"]
vocab = GeneVocab.from_file(vocab_file)
[vocab.append_token(s) for s in special_tokens if s not in vocab]

# vocab.append_token([s for s in special_tokens if s not in vocab])

logger.info("Filtering genes based on model vocab")
adata.var["id_in_vocab"] = [1 if gene in vocab else -1 for gene in genes]
    
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])

logger.info("Subsetting input data based on genes present in model vocab")
adata = adata[:, adata.var["id_in_vocab"] >= 0]

mudata.mod[par["modality"]] = adata

logger.info(f"Writing to {par['output']}")
mudata.write_h5mu(par["output"], compression=par["output_compression"])
