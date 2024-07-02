import sys
import logging
import celltypist
from mudata import read_h5mu, MuData
import anndata as ad
import scanpy as sc
import re
import numpy as np
import typing

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/tmp_TS_Blood_filtered.h5ad",
    "model": None,
    "model_save_path": "new_model.pkl",
    "only_train": False,
    "reference_obs_label": "cell_ontology_class",
    "check_expression": False,
    "feature_selection": True,
    "majority_voting": True,
    "input_obs_label": None,
    "reference_obs_label": "cell_ontology_class",
    "output_compression": "gzip"
}
meta = {}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
def setup_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger

def check_celltypist_format(indata):
    if np.abs(np.expm1(indata[0]).sum()-10000) > 1:
        return False
    return True

def main(par):
    logger = setup_logger()

    input_mudata = read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]]
        
    if par["model"]:
        model = celltypist.models.Model.load(par["model"])
    
    elif par["reference"]:
        reference_adata = ad.read_h5ad(par["reference"])
        
        reference_adata.var["gene_symbol"] = list(reference_adata.var.index)
        reference_adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_adata.var["ensemblid"]]
        
        logger.info("Detecting common vars based on ensembl ids")
        common_ens_ids = list(set(reference_adata.var.index).intersection(set(input_modality.var.index)))
        
        logger.info("  reference n_vars: %i", reference_adata.n_vars)
        logger.info("  input n_vars: %i", input_modality.n_vars)
        logger.info("  intersect n_vars: %i", len(common_ens_ids))
        assert len(common_ens_ids) >= 100, "The intersection of genes is too small."
        
        if not check_celltypist_format(input_modality.X):
            logger.warning("Input data is not in the reccommended format for CellTypist.")
        if not check_celltypist_format(reference_adata.X):
            logger.warning("Reference data is not in the reccommended format for CellTypist.")
        
        model = celltypist.train(reference_adata,
                                 labels = par["reference_obs_label"],
                                 feature_selection=par["feature_selection"],
                                 check_expression=par["check_expression"])
        
        model.write(par["model_save_path"]) if par["model_save_path"] else None
        if par["only_train"]:
            return
        
    else:
        raise ValueError("Either 'model' or 'reference' has to be provided.")
    
    predictions = celltypist.annotate(input_modality,
                                      model,
                                      majority_voting=par["majority_voting"])
    
    preds_adata = predictions.to_adata()
    input_mudata.mod[par["modality"]] = preds_adata
        
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])
    
if __name__ == '__main__':
    main(par)
