from time import sleep
import pytest
import sys
import mudata as mu
import torch
import numpy as np
from pathlib import Path
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor

## VIASH START
meta = {
    "resources_dir": "resources_test/scgpt",
    "executable": "./target/docker/scgpt/integration_pad_tokenize/integration_pad_tokenize",
    "temp_dir": "tmp",
    "config": "./target/docker/scgpt/integration_pad_tokenize/.config.vsh.yaml"
}
## VIASH END

input = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung.h5mu"
model_dir = f"{meta['resources_dir']}/scgpt/source/"
input_file = mu.read(input)

## START TEMPORARY WORKAROUND (until all scgpt modules are implemented)
# Read in data
adata = input_file.mod["rna"]

# Set tokens for integration
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]

# Make batch a category column
adata.obs["str_batch"] = adata.obs["sample"].astype(str)
batch_id_labels = adata.obs["str_batch"].astype("category").cat.codes.values
adata.obs["batch_id"] = batch_id_labels
adata.var["gene_name"] = adata.var.index.tolist()

# Load model vocab
model_dir = Path(model_dir)
vocab_file = model_dir / "vocab.json"
vocab = GeneVocab.from_file(vocab_file)
for s in special_tokens:
    if s not in vocab:
        vocab.append_token(s)

# Cross-check genes with pre-trained model
genes = adata.var["gene_name"].tolist()
adata.var["id_in_vocab"] = [
        1 if gene in vocab else -1 for gene in adata.var["gene_name"]
    ]
gene_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]

# Preprocess data
preprocessor = Preprocessor(
    use_key="X",
    filter_gene_by_counts=3,
    filter_cell_by_counts=False,
    normalize_total=10000,
    result_normed_key="X_normed",
    log1p=True,
    result_log1p_key="X_log1p",
    subset_hvg=1200,
    hvg_flavor="seurat_v3",
    binning=51,
    result_binned_key="X_binned",
    )

preprocessor(adata, batch_key="str_batch")

# copy results to mudata
input_file.mod["rna"] = adata

## END OF TEMPORARY WORKAROUND


def test_integration_pad_tokenize(run_component, tmp_path):
    output_gene_ids = tmp_path / "gene_ids.pt"
    output_values = tmp_path / "values.pt"
    output_padding_mask = tmp_path / "padding_mask.pt"
    
    input_preprocessed = f"{meta['resources_dir']}/scgpt/test_resources/Kim2020_Lung_preprocessed.h5mu"
    input_file.write(input_preprocessed)

    run_component([
        "--input", input_preprocessed,
        "--modality", "rna",
        "--output_gene_ids", output_gene_ids,
        "--output_values", output_values,
        "--output_padding_mask", output_padding_mask,
        "--pad_token", "<pad>",
        "--pad_value", "-2",
        "--input_layer", "X_binned",
        "--gene_name_layer", "gene_name",
        "--model_dir", model_dir
    ])

    adata = input_file.mod["rna"]
    gene_ids = torch.load(output_gene_ids)
    values = torch.load(output_values)
    padding_mask = torch.load(output_padding_mask)

    # check output dimensions
    ## nr of cells
    assert gene_ids.shape[0] == adata.shape[0], "gene_ids shape[0] does not match input shape[0]"
    assert values.shape[0] == adata.shape[0], "values shape[0] does not match input shape[0]"
    assert padding_mask.shape[0] == adata.shape[0], "padding_mask shape[0] does not match input shape[0]"

    ## nr of genes (subset hvg)
    assert gene_ids.shape[1] <= adata.var.shape[0] + 1, "gene_ids shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"
    assert values.shape[1] <= adata.var.shape[0] + 1, "values shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"
    assert padding_mask.shape[1] <= adata.var.shape[0] + 1, "padding_mask shape[1] is higher than adata.var.shape[0] (n_hvg + 1)"

    ## equal size of output tensors
    assert gene_ids.shape[1] == values.shape[1], "gene_ids shape[1] does not match values shape[1]"
    assert gene_ids.shape[1] == padding_mask.shape[1], "gene_ids shape[1] does not match padding_mask shape[1]"


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
