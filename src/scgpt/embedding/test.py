import pytest
import subprocess
import re
import sys
import mudata as mu
import numpy as np
from scipy.sparse import issparse
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor

## VIASH START
meta = {
    "resources_dir": "resources_test",
    "executable": "./target/docker/scgpt/integration_embedding/integration_embedding",
    "temp_dir": "tmp",
    "config": "./target/docker/scgpt/integration_embedding/.config.vsh.yaml"
}
## VIASH END

input = f"{meta['resources_dir']}/Kim2020_Lung_subset.h5mu"
model_file = f"{meta['resources_dir']}/source/best_model.pt"
vocab_file = f"{meta['resources_dir']}/source/vocab.json"
model_config_file = f"{meta['resources_dir']}/source/args.json"
input_file = mu.read(input)

## START TEMPORARY WORKAROUND DATA PREPROCESSING
#TODO: Remove this workaround once full scGPT preprocessing workflow is implemented
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
    subset_hvg=100,
    hvg_flavor="seurat_v3",
    binning=51,
    result_binned_key="X_binned",
    )

preprocessor(adata, batch_key="str_batch")

all_counts = (
    adata.layers["X_binned"].A
    if issparse(adata.layers["X_binned"])
    else adata.layers["X_binned"]
)

# Fetch gene names and look up tokens in vocab
vocab.set_default_index(vocab["<pad>"])
ntokens = len(vocab)
genes = adata.var["gene_name"].tolist()
gene_ids = np.array(vocab(genes), dtype=int)

# Fetch number of subset hvg
n_hvg = adata.var.shape[0]

# Tokenize and pad data
tokenized_data = tokenize_and_pad_batch(
    all_counts,
    gene_ids,
    max_len=n_hvg+1,
    vocab=vocab,
    pad_token=pad_token,
    pad_value=-2,
    append_cls=True,  # append <cls> token at the beginning,
    include_zero_gene=False,
    return_pt=True,
    mod_type=None,
    vocab_mod=None
    )

all_gene_ids, all_values = tokenized_data["genes"], tokenized_data["values"]
padding_mask = all_gene_ids.eq(vocab[pad_token])

adata.obsm["gene_id_tokens"] = all_gene_ids.numpy()
adata.obsm["values_tokenized"] = all_values.numpy()
adata.obsm["padding_mask"] = padding_mask.numpy()

tokenized_data = mu.MuData({'rna': adata})
tokenized_data_path = f"{meta['resources_dir']}/Kim2020_Lung_tokenized.h5mu"
tokenized_data.write_h5mu(tokenized_data_path)

## END TEMPORARY WORKAROUND DATA PREPROCESSING


def test_integration_embedding(run_component, tmp_path):

    output_embedding_file = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    run_component([
        "--input", tokenized_data_path,
        "--modality", "rna",
        "--model", model_file,
        "--model_vocab", vocab_file,
        "--model_config", model_config_file,
        "--DSBN", "True",
        "--obs_batch_label", "sample",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--obsm_padding_mask", "padding_mask",
        "--output", output_embedding_file
    ])

    # Read output file
    output_mdata = mu.read(output_embedding_file)
    output_adata = output_mdata.mod["rna"]

    # check that embedding obs is present
    assert 'X_scGPT' in output_adata.obsm.keys(), "X_scGPT is not present in anndata obsm keys"

    # check embedding size
    assert output_adata.obsm["X_scGPT"].shape[1] == 512, "Embedding size does not equal 512"

    # check embedding value range
    assert not all(np.isnan(output_adata.obsm["X_scGPT"][0])), "Embedding values are nan"
    assert all([all(i > -1) & all(i < 1) for i in output_adata.obsm["X_scGPT"]]), "Range of embedding values is outside of [-1, 1]"

    # Run embeddings without DSBN
    output_embedding_file_without_DSBN = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"
    run_component([
        "--input", tokenized_data_path,
        "--modality", "rna",
        "--model", model_file,
        "--model_vocab", vocab_file,
        "--model_config", model_config_file,
        "--DSBN", "False",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--obsm_padding_mask", "padding_mask",
        "--output", output_embedding_file_without_DSBN
    ])

    # Read output file
    output_mdata_no_DSBN = mu.read(output_embedding_file_without_DSBN)
    output_adata_no_DSBN = output_mdata_no_DSBN.mod["rna"]

    # Assert that embeddings without DSBN are different
    assert not (output_adata.obsm["X_scGPT"] == output_adata_no_DSBN.obsm["X_scGPT"]).all(), "Embeddings with and without DSBN are the same"

def test_integration_embedding_DSBN_without_batch_labels(run_component, tmp_path):
    output_embedding_file = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    args = [
        "--input", tokenized_data_path,
        "--modality", "rna",
        "--model", model_file,
        "--model_vocab", vocab_file,
        "--model_config", model_config_file,
        "--DSBN", "True",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--obsm_padding_mask", "padding_mask",
        "--output", output_embedding_file
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: When DSBN is set to True, you are required to provide batch labels \(input_obs_batch_labels\)\.",
        err.value.stdout.decode('utf-8'))

def test_integration_embedding_non_existing_keys(run_component, tmp_path):
    output_embedding_file = tmp_path / "Kim2020_Lung_subset_embedded.h5mu"

    # Test for non-existing gene names key
    args_1 = [
        "--input", tokenized_data_path,
        "--modality", "rna",
        "--model", model_file,
        "--model_vocab", vocab_file,
        "--model_config", model_config_file,
        "--DSBN", "True",
        "--obs_batch_label", "sample",
        "--var_gene_names", "dummy_gene_name_key",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--obsm_padding_mask", "padding_mask",
        "--output", output_embedding_file
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_1)
    assert re.search(
        r"KeyError: \'dummy_gene_name_key\'",
        err.value.stdout.decode('utf-8'))
    
    # Test for non-existing batch label key
    args_2 = [
        "--input", tokenized_data_path,
        "--modality", "rna",
        "--model", model_file,
        "--model_vocab", vocab_file,
        "--model_config", model_config_file,
        "--DSBN", "True",
        "--obs_batch_label", "dummy_batch_label_key",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "values_tokenized",
        "--obsm_padding_mask", "padding_mask",
        "--output", output_embedding_file
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_2)
    assert re.search(
        r"KeyError: \'dummy_batch_label_key\'",
        err.value.stdout.decode('utf-8'))

    # Test for non-existing tokenized values key
    args_3 = [
        "--input", tokenized_data_path,
        "--modality", "rna",
        "--model", model_file,
        "--model_vocab", vocab_file,
        "--model_config", model_config_file,
        "--DSBN", "True",
        "--obs_batch_label", "sample",
        "--obsm_gene_tokens", "gene_id_tokens",
        "--obsm_tokenized_values", "dummy_values_tokenized",
        "--obsm_padding_mask", "padding_mask",
        "--output", output_embedding_file
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args_3)
    assert re.search(
        r"KeyError: \'dummy_values_tokenized\'",
        err.value.stdout.decode('utf-8'))


if __name__ == '__main__':
    sys.exit(pytest.main([__file__]))
