import mudata as mu
import numpy as np
from pathlib import Path
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor


## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung.h5mu",
    "output": "resources_test/scgpt/test_resources/Kim2020_Lung_preprocessed.h5mu",
    "modality": "rna",
    "input_layer": "X",
    "ori_batch_layer_name": "sample",
    "batch_id_layer": "batch_id",
    "gene_name_layer": "gene_name",
    "normalized_total_layer": "X_normed",
    "binned_layer": "X_binned",
    "log1p_layer": "X_log1p",
    "pad_token": "<pad>",
    "filter_gene_by_counts": 3,
    "filter_cell_by_counts": -1,
    "normalize_total": 1e4,
    "n_hvg": 1200,
    "data_is_raw": False,
    "n_input_bins": 51,
    "model_dir": "resources_test/scgpt/source/",
    "output_compression": "gzip"
}
## VIASH END

# Read in data
mdata = mu.read(par["input"])
adata = mdata.mod[par["modality"]]

# Set tokens for integration
pad_token = par["pad_token"]
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Make batch a category column
adata.obs["str_batch"] = adata.obs[par["ori_batch_layer_name"]].astype(str)
batch_id_labels = adata.obs["str_batch"].astype("category").cat.codes.values
adata.obs[par["batch_id_layer"]] = batch_id_labels
adata.var[par["gene_name_layer"]] = adata.var.index.tolist()

# Load model vocab
model_dir = Path(par["model_dir"])
vocab_file = model_dir / "vocab.json"
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

# Cross-check genes with pre-trained model
genes = adata.var[par["gene_name_layer"]].tolist()
adata.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in adata.var[par["gene_name_layer"]]
    ]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]

# Set pre-processing settings
if par["filter_gene_by_counts"] == -1:
    par["filter_gene_by_counts"] = False
if par["filter_cell_by_counts"] == -1:
    par["filter_cell_by_counts"] = False
if par["normalize_total"] == -1:
    par["normalize_total"] = False

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

# copy results to mudata
mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])
