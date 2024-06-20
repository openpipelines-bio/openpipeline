import sys
import logging
import celltypist
from mudata import read_h5mu
import anndata as ad
import scanpy as sc
import re


## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/tmp_TS_Blood_filtered.h5ad",
    "model": None,
    "model_save_path": "new_model.pkl",
    "reference_obs_label": "cell_ontology_class",
    "feature_selection": True,
    "majority_voting": True
    # "prediction_mode": 
}
meta = {}
## VIASH END

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

def normalize_counts(input_adata, target_sum=1e4):
    adata = input_adata.copy()
    input_layer = adata.X
    data_for_scanpy = ad.AnnData(X=input_layer.copy())
    sc.pp.normalize_total(data_for_scanpy, target_sum=target_sum)
    adata.X = data_for_scanpy.X
    return adata

def log_transform(input_adata):
    adata = input_adata.copy()
    input_layer = adata.X
    data_for_scanpy = ad.AnnData(X=input_layer.copy())
    sc.pp.log1p(data_for_scanpy,
                base=None,
                layer=None, # use X
                copy=False) # allow overwrites in the copy that was made
    adata.X = data_for_scanpy.X
    adata.uns['log1p'] = data_for_scanpy.uns['log1p'].copy()

def main(par, meta):
    logger = setup_logger()

    input_mudata = read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]]
    
    input_modality = normalize_counts(input_modality)
    input_modality = log_transform(input_modality)

    if par["model"]:
        # TODO: implement loading of model
        pass
    
    elif par["reference"]:
        reference_adata = ad.read_h5ad(par["reference"])
        reference_adata = normalize_counts(reference_adata)
        reference_adata = log_transform(reference_adata)
        
        reference_adata.var["gene_symbol"] = list(reference_adata.var.index)
        reference_adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_adata.var["ensemblid"]]
        
        logger.info("Detecting common vars based on ensembl ids")
        common_ens_ids = list(set(reference_adata.var.index).intersection(set(input_modality.var.index)))
        
        logger.info("  reference n_vars: %i", reference_adata.n_vars)
        logger.info("  input n_vars: %i", input_modality.n_vars)
        logger.info("  intersect n_vars: %i", len(common_ens_ids))
        assert len(common_ens_ids) >= 100, "The intersection of genes is too small."
        
        model = celltypist.train(reference_adata, labels = par["reference_obs_label"], feature_selection=par["feature_selection"])
        predictions = celltypist.annotate(input_modality, model, majority_voting=par["majority_voting"])
        
    else:
        raise ValueError("Either 'model' or 'reference' has to be provided.")
