import mudata as mu
from sklearn.model_selection import train_test_split
import sys
import logging

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "test_size": 0.2,
    "val_size": None,
    "random_state": 42,
    "output_train": "train.h5mu",
    "output_val": "val.h5mu",
    "output_test": "test.h5mu",
    "compression": "gzip",
    "shuffle": True 
}
### VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()


def main():
    input_mudata = mu.read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]]
    
    n_obs = input_modality.n_obs
    train_idx, test_idx = train_test_split(range(n_obs), test_size=par["test_size"], random_state=par["random_state"], shuffle=par["shuffle"])
    
    if par["val_size"] is not None:
        val_size_relative = par["val_size"] / (1 - par["test_size"])
        train_idx, val_idx = train_test_split(train_idx, test_size=val_size_relative, random_state=par["random_state"], shuffle=par["shuffle"])
        
        train_modality = input_modality[train_idx].copy()
        val_modality = input_modality[val_idx].copy()
        test_modality = input_modality[test_idx].copy()

        train_mudata = mu.MuData({par["modality"]: train_modality})
        val_mudata = mu.MuData({par["modality"]: val_modality})
        test_mudata = mu.MuData({par["modality"]: test_modality})
        
        val_mudata.write_h5mu(par["output_val"], compression=par["compression"])
    
    else:
        train_modality = input_modality[train_idx].copy()
        test_modality = input_modality[test_idx].copy()

        train_mudata = mu.MuData({par["modality"]: train_modality})
        test_mudata = mu.MuData({par["modality"]: test_modality})
        
    train_mudata.write_h5mu(par["output_train"], compression=par["compression"])
    test_mudata.write_h5mu(par["output_test"], compression=par["compression"])


if __name__ == "__main__":
    main()
