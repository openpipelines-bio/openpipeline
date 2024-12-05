import sys
import mudata as mu
from scgpt.tokenizer.gene_tokenizer import GeneVocab

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_subset_preprocessed.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "input_var_gene_names": None,
    "output_var_filter": "id_in_vocab",
    "pad_token": "<pad>",
    "var_input": "filter_with_hvg",
    "vocab_file": "resources_test/scgpt/source/vocab.json",
    "output_compression": None,
}

meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

# Read in data
logger.info(f"Reading {par['input']}")
mudata = mu.read_h5mu(par["input"])
adata = mudata.mod[par["modality"]].copy()

pad_token = par["pad_token"]
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Fetching gene names
if not par["input_var_gene_names"]:
    genes = adata.var.index.astype(str).tolist()
elif par["input_var_gene_names"] not in adata.var.columns:
    raise ValueError(
        f"Gene name column '{par['input_var_gene_names']}' not found in .mod['{par['modality']}'].obs."
    )
else:
    genes = adata.var[par["input_var_gene_names"]].astype(str).tolist()

# Cross-check genes with pre-trained model
logger.info(f"Loading model vocab from {par['vocab_file']}")
vocab_file = par["vocab_file"]
vocab = GeneVocab.from_file(vocab_file)
[vocab.append_token(s) for s in special_tokens if s not in vocab]

if par["var_input"]:
    logger.info("Filtering genes based on model vocab and HVG")
    filter_with_hvg = adata.var[par["var_input"]].tolist()
    gene_filter_mask = [
        1 if gene in vocab and hvg else 0 for gene, hvg in zip(genes, filter_with_hvg)
    ]
    logger.info(
        f"Total number of genes after HVG present in model vocab: {str(sum(gene_filter_mask))}"
    )
else:
    logger.info("Filtering genes based on model vocab")
    gene_filter_mask = [1 if gene in vocab else 0 for gene in genes]
    logger.info(
        f"Total number of genes present in model vocab: {str(sum(gene_filter_mask))}"
    )

logger.info(f"Writing to {par['output']}")
adata.var[par["output_var_filter"]] = gene_filter_mask
adata.var[par["output_var_filter"]] = adata.var[par["output_var_filter"]].astype("bool")
mudata.mod[par["modality"]] = adata
mudata.write_h5mu(par["output"], compression=par["output_compression"])
