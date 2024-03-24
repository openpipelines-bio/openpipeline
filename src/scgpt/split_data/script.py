import sys
import os
import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.model_selection import train_test_split

## VIASH START
par = {
    "input": "resources_test/scgpt/test_resources/Kim2020_Lung.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "train_size": 0.8, 
    "test_size": 0.2,
    "val_size": 0.0,
    "shuffle": True,
    "stratify": False,
    "random_state": 42
}
## VIASH END

sys.path.append(meta["resources_dir"])
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

if not par["val_size"]:
    if par["train_size"] and par["train_size"] + par["test_size"] != 1:
        raise ValueError("The sum of all dataset partitions must be equal to 1.")
    train_obs, test_obs = train_test_split(
        adata.obs,
        test_size=par["test_size"],
        random_state=par["random_state"],
        shuffle=par["shuffle"],
        stratify=par["stratify"]
    )
    
    adata_train = adata.copy()
    adata_test = adata.copy()
    adata_train.obs = train_obs
    adata_train.obs = test_obs
    
else:
    if par["val_size"] >= 1-par["test_size"]:
        raise ValueError("Validation set size must be less than the remaining dataset size after test set partitioning.")
    if par["train_size"] and par["train_size"] + par["test_size"] + par["val_size"] != 1:
        raise ValueError("The sum of all dataset partitions must be equal to 1.")
    train_obs, test_obs = train_test_split(
        adata.obs,
        test_size=par["test_size"],
        random_state=par["random_state"],
        shuffle=par["shuffle"],
        stratify=par["stratify"]
    )
    train_obs, val_obs = train_test_split(
        train_obs,
        test_size=par["val_size"] / (1 - par["test_size"]),
        random_state=par["random_state"],
        shuffle=par["shuffle"],
        stratify=par["stratify"]
    )
    
    adata_train = adata.copy()
    adata_test = adata.copy()
    adata_val = adata.copy()
    adata_train.obs = train_obs
    adata_test.obs = test_obs
    adata_val.obs = val_obs
    
logger.info(f"Writing {par['output']}")


