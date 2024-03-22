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

