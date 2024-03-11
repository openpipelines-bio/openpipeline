import anndata as ad
import numpy as np
from pathlib import Path
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)

## VIASH START
par = {
    "input": "src/scgpt/test_resources/Kim2020_Lung.h5ad",
    "output": "src/scgpt/test_resources/Kim2020_Lung_preprocessed.h5ad",
    "input_layer": "X",
    "ori_batch_layer_name": "sample",
    "batch_id_layer": "batch_id",
    "gene_name_layer": "gene_name",
    "normalized_total_layer": "X_normed",
    "binned_layer": "X_binned",
    "log1p_layer": "X_log1p",
    "pad_token": "<pad>",
    "filter_gene_by_counts": 3,
    "filter_cell_by_counts": False,
    "normalize_total": 1e4,
    "n_hvg": 1200,
    "data_is_raw": False,
    "n_input_bins": 51,
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

# Preprocess data
preprocessor = Preprocessor(
    use_key=par["input_layer"],  # the key in adata.layers to use as raw data
    filter_gene_by_counts=par["filter_gene_by_counts"],  # step 1
    filter_cell_by_counts=par["filter_cell_by_counts"],  # step 2
    normalize_total=par["normalize_total"],  # 3. whether to normalize the raw data and to what sum
    result_normed_key=par["normalized_total_layer"],  # the key in adata.layers to store the normalized data
    log1p=par["data_is_raw"],  # 4. whether to log1p the normalized data
    result_log1p_key=par["log1p_layer"],
    subset_hvg=par["n_hvg"],  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor="seurat_v3" if par["data_is_raw"] else "cell_ranger",
    binning=par["n_input_bins"],  # 6. whether to bin the raw data and to what number of bins
    result_binned_key=par["binned_layer"],  # the key in adata.layers to store the binned data
    )

preprocessor(adata, batch_key="str_batch")

adata.write_h5ad(par["output"])
