import json
import os
import mudata as mu
from typing import Dict
import warnings
import torch
import numpy as np
from torch.utils.data import Dataset, DataLoader
from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

#TODO: check 'hyperparameter_path' (resulting in parameters) and its values
# 'embsize', 'd_hid', 'nlayers', 'nhead', 'dropout' get overridden by model_config
#TODO: subset hvg as optional

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung_preprocessed.h5mu",
    "modality": "rna",
    "input_obsm_gene_tokens": 'gene_id_tokens',
    "input_obsm_tokenized_values": 'values_tokenized',
    "model": "resources_test/scgpt/source/best_model.pt",
    # "model_config": "resources_test/scgpt/source/args.json",
    "model_vocab": "resources_test/scgpt/source/vocab.json",
    "output": "Kim2020_Lung_cell_type_annotation.h5mu",
    "gene_name_layer": "gene_name",
    "input_obs_batch_label": "sample",
    "predicted_cell_type_id": "predicted_cell_type",
    "pad_token": "<pad>",
    "dsbn": True,
    "pad_value": -2,
    "n_cls": 8,
    "n_input_bins": 51,
    "batch_size": 64,
    "output_compression": None,
    "seed": 1
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

class SeqDataset(Dataset):
    def __init__(self, data: Dict[str, torch.Tensor]):
        self.data = data

    def __len__(self):
        return self.data["gene_ids"].shape[0]

    def __getitem__(self, idx):
        return {k: v[idx] for k, v in self.data.items()}

warnings.filterwarnings('ignore')

if par["seed"]:
    set_seed(par["seed"])

logger.info(f"Setting device to {'cuda' if torch.cuda.is_available() else 'cpu'}")
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

logger.info("Reading in data")
# Read in data
mdata = mu.read(par["input"])
input_adata = mdata.mod[par["modality"]]
adata = input_adata.copy()

# Fetch batch ids for domain-specific batch normalization
if par["dsbn"] and not par["obs_batch_label"]:
    raise ValueError("When dsbn is set to True, you are required to provide batch labels (obs_batch_labels).")
elif par["dsbn"] and par["obs_batch_label"]:
    logger.info("Fetching batch id's for domain-specific batch normalization")
    batch_id_cats = adata.obs[par["obs_batch_label"]].astype("category")
    batch_id_labels = batch_id_cats.cat.codes.values
    batch_ids = batch_id_labels.tolist()
    batch_ids = np.array(batch_ids)
    num_batch_types = len(set(batch_ids))
elif not par["dsbn"] and par["obs_batch_label"]:
    logger.info("Batch labels provided but dsbn is set to False. Batch labels will be ignored and no dsbn will be performed.")

# Vocabulary configuration
special_tokens = [par["pad_token"], "<cls>", "<eoc>"]
logger.info(f"Loading model vocab from {par['model_vocab']}")
vocab_file = par["model_vocab"]
vocab = GeneVocab.from_file(vocab_file)
[vocab.append_token(s) for s in special_tokens if s not in vocab]
vocab.set_default_index(vocab[par["pad_token"]])
ntokens = len(vocab)

# Instantiate model
model = TransformerModel(
    ntokens,
    d_model=par["emb_size"], # self.encoder (GenEncoder), self.value_encoder (ContinuousValueEncoder), self.transformerencoder(TransformerEncoderLayer)
    nhead=par["n_head"], # self.transformer_encoder(TransformerEncoderLayer)
    d_hid=par["d_hid"], # self.transformer_encoder(TransformerEncoderLayer)
    nlayers=par["n_layers"], # self.transformer_encoder(TransformerEncoderLayer), self.cls_decoder
    nlayers_cls=3, # self.cls_decoder
    n_cls=par["n_cls"], # self.cls_decoder
    vocab=vocab,
    dropout=par["dropout"],  # self.transformer_encoder
    pad_token=par["pad_token"],
    pad_value=par["pad_value"],
    do_mvc=False,
    do_dab=False,
    use_batch_labels=par["dsbn"],
    num_batch_labels=num_batch_types if par["dsbn"] else None,
    domain_spec_batchnorm=par["dsbn"],
    input_emb_style="continuous",
    n_input_bins=par["n_input_bins"],
    cell_emb_style="cls",  # hard-coded
    use_fast_transformer=False,   #TODO: parametrize when GPU is available
    fast_transformer_backend="flash",  #TODO: parametrize when GPU is available
    pre_norm=False,  #TODO: parametrize when GPU is available
)

model_file = par["model"]
try:
    logger.info(f"Loading all model params from {model_file}")
    model.load_state_dict(torch.load(model_file, map_location=device))
except:
    logger.info("only load params that are in the model and match the size")
    model_dict = model.state_dict()
    pretrained_dict = torch.load(model_file, map_location=device)
    pretrained_dict = {
        k: v
        for k, v in pretrained_dict.items()
        if k in model_dict and v.shape == model_dict[k].shape
    }
    for k, v in pretrained_dict.items():
        logger.info(f"Loading params {k} with shape {v.shape}")
    model_dict.update(pretrained_dict)
    model.load_state_dict(model_dict)

model.to(device)

input_gene_ids = adata.obsm[par["obsm_gene_tokens"]]
input_values = adata.obsm[par["obsm_tokenized_values"]]

data_pt = {
    "gene_ids": input_gene_ids,
    "values": input_values,
    "batch_labels": torch.from_numpy(batch_ids).long(),
}

data_loader = DataLoader(
    dataset=SeqDataset(data_pt),
    batch_size=par["batch_size"],
    num_workers=min(len(os.sched_getaffinity(0)), par["batch_size"] // 2),
    pin_memory=True,
)

logger.info("Predicting cell type labels")
model.eval()
predictions = []
with torch.no_grad():
    for batch_data in data_loader:
        input_gene_ids = batch_data["gene_ids"].to(device)
        input_values = batch_data["values"].to(device)
        batch_labels = batch_data["batch_labels"].to(device)

        src_key_padding_mask = input_gene_ids.eq(vocab[par["pad_token"]])
        with torch.cuda.amp.autocast(enabled=False): # parametrize
            output_dict = model(
                input_gene_ids,
                input_values,
                src_key_padding_mask=src_key_padding_mask,
                batch_labels=batch_labels if par["dsbn"] else None,
                CLS=True  # Return celltype classification objective output
            )
            output_values = output_dict["cls_output"]
        preds = output_values.argmax(1).cpu().numpy()
        predictions.append(preds)

predictions = np.concatenate(predictions, axis=0)
adata.obs[par["obs_predicted_cell_type"]] = predictions

# Write output
logger.info("Writing output data")
mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])
