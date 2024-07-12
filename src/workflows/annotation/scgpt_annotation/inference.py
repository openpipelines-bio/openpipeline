import copy
import gc
import yaml
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback
from typing import List, Tuple, Dict, Union, Optional
import warnings
import pandas as pd
import pickle
import torch
from anndata import AnnData
import scanpy as sc
import scvi
import seaborn as sns
import numpy as np
import wandb
from scipy.sparse import issparse
import matplotlib.pyplot as plt
from torch import nn
from torch.nn import functional as F
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from torchtext.vocab import Vocab
from torchtext._torchtext import (
    Vocab as VocabPybind,
)
from sklearn.metrics import confusion_matrix
import argparse
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from torch.nn.functional import softmax

parser = argparse.ArgumentParser()
parser.add_argument("hyperparameter_path", type=str, help="The path to the JSON file containing the hyperparameters")
args = parser.parse_args()

with open(args.hyperparameter_path, 'r') as file:
    parameters = yaml.safe_load(file)
run_params = copy.deepcopy(parameters)


test_dataset = parameters.pop('test_dataset')
scgpt_path = parameters.pop('scgpt_path')
target_column_test = parameters.pop('target_column_test')
apply_mapping_dic = parameters.pop('apply_mapping_dic')


sys.path.insert(0, scgpt_path)
import scgpt as scg
from scgpt.model import TransformerModel, AdversarialDiscriminator
from scgpt.tokenizer import tokenize_and_pad_batch, random_mask_value
from scgpt.loss import (
    masked_mse_loss,
    masked_relative_error,
    criterion_neg_log_bernoulli,
)

from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.preprocess import Preprocessor
from scgpt import SubsetsBatchSampler
from scgpt.utils import set_seed, category_str2int, eval_scib_metrics

class SeqDataset(Dataset):
    def __init__(self, data: Dict[str, torch.Tensor]):
        self.data = data

    def __len__(self):
        return self.data["gene_ids"].shape[0]

    def __getitem__(self, idx):
        return {k: v[idx] for k, v in self.data.items()}

sc.set_figure_params(figsize=(6, 6))
os.environ["KMP_WARNINGS"] = "off"
warnings.filterwarnings('ignore')

run = wandb.init(
    config=parameters,
    project="scGPT",
    reinit=True,
    settings=wandb.Settings(start_method="fork"),
)
config = wandb.config
print(config)

set_seed(config.seed)

# settings for input and preprocessing
pad_token = "<pad>"
special_tokens = [pad_token, "<cls>", "<eoc>"]
mask_ratio = config.mask_ratio
mask_value = "auto"  # for masked values, now it should always be auto

include_zero_gene = config.include_zero_gene  # if True, include zero genes among hvgs in the training
max_seq_len = 3001
n_bins = config.n_bins

# input/output representation
input_style = "binned"  # "normed_raw", "log1p", or "binned"
output_style = "binned"  # "normed_raw", "log1p", or "binned"

# settings for training
MLM = False  # whether to use masked language modeling, currently it is always on.
CLS = True  # celltype classification objective
ADV = False  # Adversarial training for batch correction
CCE = False  # Contrastive cell embedding objective
MVC = config.MVC  # Masked value prediction for cell embedding
ECS = config.ecs_thres > 0  # Elastic cell similarity objective
DAB = False  # Domain adaptation by reverse backpropagation, set to 2 for separate optimizer
INPUT_BATCH_LABELS = False  # TODO: have these help MLM and MVC, while not to classifier
input_emb_style = "continuous"  # "category" or "continuous" or "scaling"
cell_emb_style = "cls"  # "avg-pool" or "w-pool" or "cls"
adv_E_delay_epochs = 0  # delay adversarial training on encoder for a few epochs
adv_D_delay_epochs = 0
mvc_decoder_style = "inner product"
ecs_threshold = config.ecs_thres
dab_weight = config.dab_weight

explicit_zero_prob = MLM and include_zero_gene  # whether explicit bernoulli for zeros
do_sample_in_train = False and explicit_zero_prob  # sample the bernoulli in training

per_seq_batch_sample = False

# settings for optimizer
lr = config.lr  # TODO: test learning rate ratio between two tasks
lr_ADV = 1e-3  # learning rate for discriminator, used when ADV is True
batch_size = config.batch_size
eval_batch_size = config.batch_size
epochs = config.epochs
schedule_interval = 1

# settings for the model
fast_transformer = config.fast_transformer
fast_transformer_backend = "flash"  # "linear" or "flash"
embsize = config.layer_size  # embedding dimension
d_hid = config.layer_size  # dimension of the feedforward network in TransformerEncoder
nlayers = config.nlayers  # number of TransformerEncoderLayer in TransformerEncoder
nhead = config.nhead  # number of heads in nn.MultiheadAttention
dropout = config.dropout  # dropout probability

# logging
log_interval = 100  # iterations
save_eval_interval = config.save_eval_interval  # epochs
do_eval_scib_metrics = True

assert input_style in ["normed_raw", "log1p", "binned"]
assert output_style in ["normed_raw", "log1p", "binned"]
assert input_emb_style in ["category", "continuous", "scaling"]
if input_style == "binned":
    if input_emb_style == "scaling":
        raise ValueError("input_emb_style `scaling` is not supported for binned input.")
elif input_style == "log1p" or input_style == "normed_raw":
    if input_emb_style == "category":
        raise ValueError(
            "input_emb_style `category` is not supported for log1p or normed_raw input."
        )

if input_emb_style == "category":
    mask_value = n_bins + 1
    pad_value = n_bins  # for padding gene expr values
    n_input_bins = n_bins + 2
else:
    mask_value = -1
    pad_value = -2
    n_input_bins = n_bins

if ADV and DAB:
    raise ValueError("ADV and DAB cannot be both True.")
DAB_separate_optim = True if DAB > 1 else False

dataset_name = config.dataset_name
save_dir = Path(f"./save/test/{dataset_name}_{time.strftime('%b%d-%H-%M')}/")
save_dir.mkdir(parents=True, exist_ok=True)
print(f"save to {save_dir}")
logger = scg.logger
scg.utils.add_file_handler(logger, save_dir / "run.log")


saved_model_dir_path = Path(config.load_model)
saved_model_file = saved_model_dir_path / "best_model.pt"
# Load the checkpoint
saved_checkpoint = torch.load(saved_model_file)
#Extract the mapping dictionary 
mapping_dic = saved_checkpoint['mapping_dic']
print("mapping_dic:",mapping_dic)
# Extract the id2type dictionary
id2type = saved_checkpoint['id_to_class']
print("id2type:",id2type)

#read test dataset
adata_test = sc.read_h5ad(test_dataset) # test set

#apply mapping dictionary
if apply_mapping_dic:
    adata_test.obs[target_column_test] = adata_test.obs[target_column_test].map(mapping_dic)

print(adata_test.obs[target_column_test].value_counts())

adata_test.obs["str_batch"] = "1"
adata_test.obs["celltype"] = adata_test.obs[target_column_test].astype("category")
adata_test.obs["batch_id"]  = adata_test.obs["str_batch"] = "1"

data_is_raw = False
filter_gene_by_counts = False
adata_test_raw = adata_test.copy()

batch_id_labels = adata_test.obs['str_batch'].astype("category").cat.codes.values
adata_test.obs["batch_id"] = batch_id_labels
celltypes = adata_test.obs["celltype"].unique()
num_types = len(celltypes)
print('num_types:',num_types)


type2id = {v: k for k, v in id2type.items()}

# apply the mapping to test dataset
default_id = max(type2id.values()) + 1  # default ID for unknown/new cell types
adata_test.obs["celltype_id"] = adata_test.obs[target_column_test].apply(lambda x: type2id.get(x, default_id))

adata_test.var["gene_name"] = adata_test.var.index.tolist()
batch_ids = adata_test.obs["batch_id"].tolist()
num_batch_types = len(set(batch_ids))

if config.load_model is not None:
    model_dir = Path(config.load_model)
    model_config_file = model_dir / "args.json"
    model_file = model_dir / "best_model.pt"
    vocab_file = model_dir / "vocab.json"

    vocab = GeneVocab.from_file(vocab_file)
    shutil.copy(vocab_file, save_dir / "vocab.json")
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
            
    # model
    with open(model_config_file, "r") as f:
        model_configs = json.load(f)
    logger.info(
        f"Resume model from {model_file}, the model args will override the "
        f"config {model_config_file}."
    )
    embsize = model_configs["embsize"]
    nhead = model_configs["nheads"]
    d_hid = model_configs["d_hid"]
    nlayers = model_configs["nlayers"]
    n_layers_cls = model_configs["n_layers_cls"]

preprocessor = Preprocessor(
    use_key="X",  # the key in adata.layers to use as raw data
    filter_gene_by_counts=filter_gene_by_counts,  # step 1
    filter_cell_by_counts=False,  # step 2
    normalize_total=1e4,  # 3. whether to normalize the raw data and to what sum
    result_normed_key="X_normed",  # the key in adata.layers to store the normalized data
    log1p=data_is_raw,  # 4. whether to log1p the normalized data
    result_log1p_key="X_log1p",
    subset_hvg=False,  # 5. whether to subset the raw data to highly variable genes
    hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
    binning=n_bins,  # 6. whether to bin the raw data and to what number of bins
    result_binned_key="X_binned",  # the key in adata.layers to store the binned data
)

preprocessor(adata_test, batch_key=None)

input_layer_key = {  # the values of this map coorespond to the keys in preprocessing
    "normed_raw": "X_normed",
    "log1p": "X_normed",
    "binned": "X_binned",
}[input_style]

genes = adata_test.var["gene_name"].tolist()

if config.load_model is None:
    vocab = Vocab(
        VocabPybind(genes + special_tokens, None)
    )  # bidirectional lookup [gene <-> int]
vocab.set_default_index(vocab["<pad>"])
gene_ids = np.array(vocab(genes), dtype=int)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
ntokens = len(vocab)  # size of vocabulary
model = TransformerModel(
    ntokens,
    embsize,
    nhead,
    d_hid,
    nlayers,
    nlayers_cls=3,
    n_cls=num_types if CLS else 1,
    vocab=vocab,
    dropout=dropout,
    pad_token=pad_token,
    pad_value=pad_value,
    do_mvc=MVC,
    do_dab=DAB,
    use_batch_labels=INPUT_BATCH_LABELS,
    num_batch_labels=num_batch_types,
    domain_spec_batchnorm=config.DSBN,
    input_emb_style=input_emb_style,
    n_input_bins=n_input_bins,
    cell_emb_style=cell_emb_style,
    mvc_decoder_style=mvc_decoder_style,
    ecs_threshold=ecs_threshold,
    explicit_zero_prob=explicit_zero_prob,
    use_fast_transformer=fast_transformer,
    fast_transformer_backend=fast_transformer_backend,
    pre_norm=config.pre_norm,
)
if config.load_model is not None:
    try:
        checkpoint = torch.load(model_file)
        model.load_state_dict(checkpoint['model_state_dict'])
        #model.load_state_dict(torch.load(model_file))
        logger.info(f"Loading all model params from {model_file}")
    except:
        # only load params that are in the model and match the size
        model_dict = model.state_dict()
        pretrained_dict = torch.load(model_file)
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
wandb.watch(model)

if ADV:
    discriminator = AdversarialDiscriminator(
        d_model=embsize,
        n_cls=num_batch_types,
    ).to(device)
criterion = masked_mse_loss
criterion_cls = nn.CrossEntropyLoss()
criterion_dab = nn.CrossEntropyLoss()
optimizer = torch.optim.Adam(
    model.parameters(), lr=lr, eps=1e-4 if config.amp else 1e-8
)
scheduler = torch.optim.lr_scheduler.StepLR(
    optimizer, schedule_interval, gamma=config.schedule_ratio
)
if DAB_separate_optim:
    optimizer_dab = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler_dab = torch.optim.lr_scheduler.StepLR(
        optimizer_dab, schedule_interval, gamma=config.schedule_ratio
    )
if ADV:
    criterion_adv = nn.CrossEntropyLoss()  # consider using label smoothing
    optimizer_E = torch.optim.Adam(model.parameters(), lr=lr_ADV)
    scheduler_E = torch.optim.lr_scheduler.StepLR(
        optimizer_E, schedule_interval, gamma=config.schedule_ratio
    )
    optimizer_D = torch.optim.Adam(discriminator.parameters(), lr=lr_ADV)
    scheduler_D = torch.optim.lr_scheduler.StepLR(
        optimizer_D, schedule_interval, gamma=config.schedule_ratio
    )

scaler = torch.cuda.amp.GradScaler(enabled=config.amp)
def define_wandb_metrcis():
    wandb.define_metric("valid/mse", summary="min", step_metric="epoch")
    wandb.define_metric("valid/mre", summary="min", step_metric="epoch")
    wandb.define_metric("valid/dab", summary="min", step_metric="epoch")
    wandb.define_metric("valid/sum_mse_dab", summary="min", step_metric="epoch")
    wandb.define_metric("test/avg_bio", summary="max")

def evaluate(model: nn.Module, loader: DataLoader, return_raw: bool = False) -> float:
    """
    Evaluate the model on the evaluation data.
    """
    model.eval()
    total_loss = 0.0
    total_error = 0.0
    total_dab = 0.0
    total_num = 0
    predictions = []
    confidences = []  # Store confidence scores
    probabilities = []  # Store probabilities for all classes
    with torch.no_grad():
        for batch_data in loader:
            input_gene_ids = batch_data["gene_ids"].to(device)
            input_values = batch_data["values"].to(device)
            target_values = batch_data["target_values"].to(device)
            batch_labels = batch_data["batch_labels"].to(device)
            celltype_labels = batch_data["celltype_labels"].to(device)

            src_key_padding_mask = input_gene_ids.eq(vocab[pad_token])
            with torch.cuda.amp.autocast(enabled=config.amp):
                output_dict = model(
                    input_gene_ids,
                    input_values,
                    src_key_padding_mask=src_key_padding_mask,
                    batch_labels=batch_labels if INPUT_BATCH_LABELS or config.DSBN else None,
                    CLS=CLS,  # evaluation does not need CLS or CCE
                    CCE=False,
                    MVC=False,
                    ECS=False,
                    do_sample=do_sample_in_train,
                    #generative_training = False,
                )
                output_values = output_dict["cls_output"]
                loss = criterion_cls(output_values, celltype_labels)

                if DAB:
                    loss_dab = criterion_dab(output_dict["dab_output"], batch_labels)
            
            # Compute softmax to get probabilities
            probs = F.softmax(output_values, dim=1)
            probabilities.append(probs.cpu().numpy())
            
            total_loss += loss.item() * len(input_gene_ids)
            accuracy = (output_values.argmax(1) == celltype_labels).sum().item()
            total_error += (1 - accuracy / len(input_gene_ids)) * len(input_gene_ids)
            total_dab += loss_dab.item() * len(input_gene_ids) if DAB else 0.0
            total_num += len(input_gene_ids)
            preds = output_values.argmax(1).cpu().numpy()
            predictions.append(preds)
            # Get confidence scores
            confs = softmax(output_values, dim=1).max(1)[0].cpu().numpy()  # Get max probability as confidence
            confidences.append(confs)

    wandb.log(
        {
            "valid/mse": total_loss / total_num,
            "valid/err": total_error / total_num,
            "valid/dab": total_dab / total_num,
            "valid/sum_mse_dab": (total_loss + dab_weight * total_dab) / total_num,
        },
    )
    if return_raw:
        # Concatenate lists of numpy arrays along the first axis
        predictions_concat = np.concatenate(predictions, axis=0)
        confidences_concat = np.concatenate(confidences, axis=0)
        probabilities_concat = np.concatenate(probabilities, axis=0)
        return predictions_concat, confidences_concat, probabilities_concat
        
    return total_loss / total_num, total_error / total_num



def test(model: nn.Module, adata: DataLoader) -> float:
    all_counts = (
        adata.layers[input_layer_key].A
        if issparse(adata.layers[input_layer_key])
        else adata.layers[input_layer_key]
    )

    celltypes_labels = adata.obs["celltype_id"].tolist()  # make sure count from 0
    celltypes_labels = np.array(celltypes_labels)

    batch_ids = adata.obs["batch_id"].tolist()
    batch_ids = np.array(batch_ids)

    tokenized_test = tokenize_and_pad_batch(
        all_counts,
        gene_ids,
        max_len=max_seq_len,
        vocab=vocab,
        pad_token=pad_token,
        pad_value=pad_value,
        append_cls=True,  # append <cls> token at the beginning
        include_zero_gene=include_zero_gene,
    )

    input_values_test = random_mask_value(
        tokenized_test["values"],
        mask_ratio=mask_ratio,
        mask_value=mask_value,
        pad_value=pad_value,
    )

    test_data_pt = {
        "gene_ids": tokenized_test["genes"],
        "values": input_values_test,
        "target_values": tokenized_test["values"],
        "batch_labels": torch.from_numpy(batch_ids).long(),
        "celltype_labels": torch.from_numpy(celltypes_labels).long(),
    }

    test_loader = DataLoader(
        dataset=SeqDataset(test_data_pt),
        batch_size=eval_batch_size,
        shuffle=False,
        drop_last=False,
        num_workers=min(len(os.sched_getaffinity(0)), eval_batch_size // 2),
        pin_memory=True,
    )

    model.eval()
    predictions,confidence,probabilities = evaluate(
        model,
        loader=test_loader,
        return_raw=True,
    )
        # compute accuracy, precision, recall, f1
    top_probabilities = probabilities.max(axis=1)
    accuracy = accuracy_score(celltypes_labels, predictions)
    precision = precision_score(celltypes_labels, predictions, average="macro")
    recall = recall_score(celltypes_labels, predictions, average="macro")
    macro_f1 = f1_score(celltypes_labels, predictions, average="macro")

    logger.info(
        f"Accuracy: {accuracy:.3f}, Precision: {precision:.3f}, Recall: {recall:.3f}, "
        f"Macro F1: {macro_f1:.3f}"
    )

    results = {
        "test/accuracy": accuracy,
        "test/precision": precision,
        "test/recall": recall,
        "test/macro_f1": macro_f1,
    }

    return confidence, predictions, celltypes_labels, results, probabilities, top_probabilities

confidences,predictions, labels, results,probabilities,top_probabilities = test(model, adata_test)

id2type[default_id] = "Unknown"

adata_test_raw.obs["predictions"] = [id2type[p] for p in predictions]
adata_test_raw.obs["confidence_scores"] = confidences  # Store confidence scores

celltypes = sorted(list(set(adata_test_raw.obs['celltype']) | set(adata_test_raw.obs['predictions'])))

# Create a sorted palette
palette_ = plt.rcParams["axes.prop_cycle"].by_key()["color"] * 3  # Ensure enough colors
palette_ = {c: palette_[i] for i, c in enumerate(celltypes)}

# Ensure the categories in the adata object are ordered
adata_test_raw.obs['celltype'] = adata_test_raw.obs['celltype'].astype('category')
adata_test_raw.obs['celltype'].cat.set_categories(celltypes, ordered=True, inplace=True)

adata_test_raw.obs['predictions'] = adata_test_raw.obs['predictions'].astype('category')
adata_test_raw.obs['predictions'].cat.set_categories(celltypes, ordered=True, inplace=True)

# Create the figure with subplots
fig, axs = plt.subplots(2, 1, figsize=(12, 12), dpi=300)  # Adjust the figsize if necessary

# Plot the test UMAP and predictions UMAP
sc.pl.umap(adata_test_raw, color="celltype", palette=palette_, show=False, ax=axs[0], title="Test UMAP")
sc.pl.umap(adata_test_raw, color="predictions", palette=palette_, show=False, ax=axs[1], title="Predictions UMAP")

plt.tight_layout()
plt.savefig(save_dir / "umap.png", dpi=300)

'''
Correct vs Incorrect UMAP
'''

# Creating a new column for correct vs incorrect predictions
adata_test_raw.obs['correct'] = (adata_test_raw.obs['celltype'] == adata_test_raw.obs['predictions'])
adata_test_raw.obs['correct'] = adata_test_raw.obs['correct'].map({True: 'Correct', False: 'Incorrect'})

# Create a custom color palette for correct and incorrect predictions
palette_correct_incorrect = {'Correct': 'grey', 'Incorrect': 'red'}

# Create the figure with subplots
fig, axs = plt.subplots(2, 1, figsize=(12, 12), dpi=300)  # Adjust the figsize if necessary

# Plot the UMAP for test data (cell types)
sc.pl.umap(adata_test_raw, color='celltype', palette=palette_, show=False, ax=axs[0], title='Test UMAP')

# Plot the UMAP for correct vs incorrect predictions
sc.pl.umap(adata_test_raw, color='correct', palette=palette_correct_incorrect, show=False, ax=axs[1], title='Correct vs Incorrect Predictions UMAP')

# Ensure a tight layout and save the figure
plt.tight_layout()
plt.savefig(save_dir / 'umap_correct_incorrect.png', dpi=300)

# Saving results and additional visualizations
save_dict = {
    "predictions": predictions,
    "confidences":confidences,
    "labels": labels,
    "results": results,
    "id_maps": id2type,
    "probabilites":probabilities,
    "top_probabilites":top_probabilities
}

with open(save_dir / "results.pkl", "wb") as f:
    pickle.dump(save_dict, f)

# Log UMAP with WandB
results["test/cell_umap"] = wandb.Image(str(save_dir / "umap.png"), caption=f"predictions macro f1 {results['test/macro_f1']:.3f}")
wandb.log(results)

# Normalize and plot the confusion matrix
labels_mapped = [id2type[l] for l in labels]  # Map numeric labels to cell type names
predictions_mapped = adata_test_raw.obs["predictions"].tolist()  # Use mapped predictions

# Adjust cell types list based on available data for true labels
celltypes_true = sorted(set(labels_mapped))

# Adjust cell types list based on actual predictions made
celltypes_pred = sorted(set(predictions_mapped))

# Compute confusion matrix without normalization to filter predicted labels afterwards
cm = confusion_matrix(labels_mapped, predictions_mapped, labels=celltypes_true,normalize='true')

# Filtering out columns from the confusion matrix that correspond to predicted labels without any predictions
cm_filtered = cm[:, [celltypes_true.index(pred) for pred in celltypes_pred if pred in celltypes_true]]

# Ensure the DataFrame matches the filtered confusion matrix shape
cm_df = pd.DataFrame(cm_filtered, index=celltypes_true, columns=[pred for pred in celltypes_pred if pred in celltypes_true])

plt.figure(figsize=(18, 18))
sns.heatmap(cm_df, annot=True, fmt=".2f", cmap="Blues", xticklabels=True, yticklabels=True)
plt.ylabel('True Label')
plt.xlabel('Predicted Label')
plt.savefig(save_dir / "confusion_matrix.png", dpi=300)

# Processing confidence scores for specific cell types and overall distribution
confidence_scores = adata_test_raw.obs["confidence_scores"]

# Summary of confidence scores
print("Confidence Scores Summary:")
print(f"Min: {confidence_scores.min()}, Max: {confidence_scores.max()}, Mean: {confidence_scores.mean()}, Median: {np.median(confidence_scores)}")

plt.figure(figsize=(18, 18))
sns.histplot(confidence_scores, kde=True, color='blue', bins=30)
plt.title('Distribution of Confidence Scores')
plt.xlabel('Confidence Score')
plt.ylabel('Frequency')
plt.savefig(save_dir / "confidence.png", dpi=300)

# Log confusion matrix with WandB
results["test/confusion_matrix"] = wandb.Image(str(save_dir / "confusion_matrix.png"), caption="confusion matrix")
wandb.log(results)

# Saving test run parameters

with open(save_dir / "test_params.json", "w") as file:
    json.dump(run_params, file)