import json
import os
import mudata as mu
from typing import Dict
import warnings
import torch
import numpy as np
from torch.nn import functional
from torch.utils.data import Dataset, DataLoader
from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

## VIASH START
par = {
  'input': r'resources_test/scgpt/test_resources/Kim2020_Lung_subset_tokenized.h5mu',
  'modality': r'rna',
  'model': r'resources_test/scgpt/source/best_model.pt',
  'model_config': r'resources_test/scgpt/source/args.json',
  'model_vocab': r'resources_test/scgpt/source/vocab.json',
  'obs_batch_label': r'sample',
  'obsm_gene_tokens': r'gene_id_tokens',
  'obsm_tokenized_values': r'values_tokenized',
  'output': r'output.h5mu',
  'output_compression': None,
  'obs_predicted_cell_class': r'predicted_cell_class',
  'obs_predicted_cell_label': r'predicted_cell_label',
  'dsbn': True,
  'seed': 0,
  'pad_token': "<pad>",
  'pad_value': -2,
  'n_input_bins': 51,
  'batch_size': 64,
  'finetuned_checkpoints_key': 'mapping_dic',
  'label_mapper_key': 'id_to_class'
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

# Setting seed
if par["seed"]:
    set_seed(par["seed"])

# Setting device
logger.info(f"Setting device to {'cuda' if torch.cuda.is_available() else 'cpu'}")
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# Read in data
logger.info("Reading in data")
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
elif not par["dsbn"]:
    # forward pass requires a tensor as input
    batch_ids = np.zeros(adata.shape[0])

# Vocabulary configuration
logger.info("Loading model vocabulary")
special_tokens = [par["pad_token"], "<cls>", "<eoc>"]
logger.info(f"Loading model vocab from {par['model_vocab']}")
vocab_file = par["model_vocab"]
vocab = GeneVocab.from_file(vocab_file)
[vocab.append_token(s) for s in special_tokens if s not in vocab]
vocab.set_default_index(vocab[par["pad_token"]])
ntokens = len(vocab)

# Model configuration
logger.info("Loading model and configurations")
model_config_file = par["model_config"]
with open(model_config_file, "r") as f:
    model_configs = json.load(f)
embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]

# Ensure the provided model has the correct architecture
logger.info("Loading model")
model_file = par["model"]
model_dict = torch.load(model_file, map_location=device)
for k, v in {
        "--finetuned_checkpoints_key": par["finetuned_checkpoints_key"],
        "--label_mapper_key": par["label_mapper_key"],
        }.items():
    if v not in model_dict.keys():
        raise KeyError(f"The key '{v}' provided for '{k}' could not be found in the provided --model file. The finetuned model file for cell type annotation requires valid keys for the checkpoints and the label mapper.")
pretrained_dict = model_dict[par["finetuned_checkpoints_key"]]

# Label mapper configuration
logger.info("Loading label mapper")
label_mapper = model_dict[par["label_mapper_key"]]
cell_type_mapper = {int(k): v for k, v in label_mapper.items()}
n_cls = len(cell_type_mapper)

# Model instatiation
logger.info("Instantiating model")
model = TransformerModel(
    ntokens,
    d_model=embsize,  # self.encoder (GenEncoder), self.value_encoder (ContinuousValueEncoder), self.transformerencoder(TransformerEncoderLayer)
    nhead=nhead,  # self.transformer_encoder(TransformerEncoderLayer)
    d_hid=d_hid,  # self.transformer_encoder(TransformerEncoderLayer)
    nlayers=nlayers,  # self.transformer_encoder(TransformerEncoderLayer), self.cls_decoder
    nlayers_cls=3,  # self.cls_decoder
    n_cls=n_cls,  # self.cls_decoder
    vocab=vocab,
    dropout=0.2,  # self.transformer_encoder
    pad_token=par["pad_token"],
    pad_value=par["pad_value"],
    do_mvc=False,
    do_dab=False,
    use_batch_labels=par["dsbn"],
    num_batch_labels=num_batch_types if par["dsbn"] else None,
    domain_spec_batchnorm=par["dsbn"],
    input_emb_style="continuous",
    n_input_bins=par["n_input_bins"],
    cell_emb_style="cls",  # required for cell-type annotation
    use_fast_transformer=False,   #TODO: parametrize when GPU is available
    fast_transformer_backend="flash",  #TODO: parametrize when GPU is available
    pre_norm=False,  #TODO: parametrize when GPU is available
)


# Load model params
logger.info(f"Loading model params from {model_file}")
try:
    model.load_state_dict(pretrained_dict)
except RuntimeError:
    logger.info("only load params that are in the model and match the size")
    model_dict = model.state_dict()
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

# Load tokenized gene data
logger.info("Loading data for inference")
for k, v in {
        "--obsm_gene_tokens": par["obsm_gene_tokens"],
        "--obsm_tokenized_values": par["obsm_tokenized_values"],
        }.items():
    if v not in adata.obsm.keys():
        raise KeyError(f"The parameter '{v}' provided for '{k}' could not be found in adata.obsm")

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
    num_workers=min(os.cpu_count(), par["batch_size"] // 2),
    pin_memory=True,
)

# Inference
logger.info("Predicting cell type classes")
model.eval()
predictions = []
probabilities = []
confidences = []
with torch.no_grad():
    for batch_data in data_loader:
        input_gene_ids = batch_data["gene_ids"].to(device)
        input_values = batch_data["values"].to(device)
        batch_labels = batch_data["batch_labels"].to(device)

        src_key_padding_mask = input_gene_ids.eq(vocab[par["pad_token"]])
        with torch.cuda.amp.autocast(enabled=False):
            output_dict = model(
                input_gene_ids,
                input_values,
                src_key_padding_mask=src_key_padding_mask,
                batch_labels=batch_labels if par["dsbn"] else None,
                CLS=True,  # Return celltype classification objective output
                CCE=False,
                MVC=False,
                ECS=False,
            )
            output_values = output_dict["cls_output"]

        preds = output_values.argmax(1).cpu().numpy()
        predictions.append(preds)

        probs = functional.softmax(output_values, dim=1).max(1)[0]
        probabilities.append(probs.cpu().numpy())

predictions = np.concatenate(predictions, axis=0)
probabilities = np.concatenate(probabilities, axis=0)

# Assign cell type labels to predicted classes
logger.info("Assigning cell type predictions and probabilities")
adata.obs["scgpt_class_pred"] = predictions
adata.obs[par["output_obs_predictions"]] = adata.obs["scgpt_class_pred"].map(lambda x: cell_type_mapper[x])
adata.obs[par["output_obs_probability"]] = probabilities

# Write output
logger.info("Writing output data")
mdata.mod[par["modality"]] = adata
mdata.write(par["output"], compression=par["output_compression"])