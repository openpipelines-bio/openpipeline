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
    "reference_obs_label": "cell_ontology_class",
    "feature_selection": True,
    "majority_voting": True,
    "input_layer": None,
    "input_obs_label": None,
    "input_obs_batch": None,
    "input_var_subset": None,
    "reference_layer": None,
    "reference_obs_label": "cell_ontology_class",
    "reference_obs_batch": "donor_assay",
    "output_compression": "gzip"
    # "prediction_mode": 
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

# Helper functions
def get_X(adata: ad.AnnData, layer: typing.Optional[str], var_index: typing.Optional[str]):
    """Fetch the counts data from X or a layer. Subset columns by var_index if so desired."""
    if var_index:
        adata = adata[:, var_index]
    if layer:
        return adata.layers[layer]
    else:
        return adata.X
def get_obs(adata: ad.AnnData, obs_par_names):
    """Subset the obs dataframe to just the columns defined by the obs_label and obs_batch."""
    obs_columns = [par[x] for x in obs_par_names if par[x]]
    return adata.obs[obs_columns]
def get_var(adata: ad.AnnData, var_index: list[str]):
    """Fetch the var dataframe. Subset rows by var_index if so desired."""
    return adata.var.loc[var_index]

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
        
        input_modality = ad.AnnData(
            X = get_X(input_modality, par["input_layer"], common_ens_ids),
            obs = get_obs(input_modality, ["input_obs_label", "input_obs_batch"]),
            var = get_var(input_modality, common_ens_ids)
        )
        reference_adata = ad.AnnData(
            X = get_X(reference_adata, par["reference_layer"], common_ens_ids),
            obs = get_obs(reference_adata, ["reference_obs_label", "reference_obs_batch"]),
            var = get_var(reference_adata, common_ens_ids)
        )
        
        if not check_celltypist_format(input_modality.X):
            logger.warning("Input data is not in the reccommended format for CellTypist.")
        if not check_celltypist_format(reference_adata.X):
            logger.warning("Reference data is not in the reccommended format for CellTypist.")
        
        model = celltypist.train(reference_adata,
                                 labels = par["reference_obs_label"],
                                 feature_selection=par["feature_selection"],
                                 check_expression=False)
        
    else:
        raise ValueError("Either 'model' or 'reference' has to be provided.")
    
    predictions = celltypist.annotate(input_modality, model, majority_voting=par["majority_voting"])
    preds_adata = predictions.to_adata()
    input_mudata.mod[par["modality"]] = preds_adata
        
    model.write(par["model_save_path"]) if par["model_save_path"] else None
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])
    
if __name__ == '__main__':
    main(par)
