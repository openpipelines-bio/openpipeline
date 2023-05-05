import logging
import mudata
import numpy as np
from pandas import DataFrame
import scvi
from torch.cuda import is_available as cuda_is_available
try:
    from torch.backends.mps import is_available as mps_is_available
except ModuleNotFoundError:
    # Older pytorch versions
    # MacOS GPUs
    def mps_is_available():
        return False
    
import os

### VIASH START
par = {
    "query": "resources_test/scvi_tools/totalvi_pbmc_query_test.h5mu",
    "reference": "resources_test/scvi_tools/totalvi_pbmc_ref_test.h5mu",
    "query_modality": "rna",
    "query_proteins_key": None,
    "reference_modality": "rna",
    "reference_proteins_key": None,
    "force_retrain": False,
    "input_layer": None,
    "obs_batch": "sample_id",
    "var_input": None,
    "output": "foo.h5mu",
    "obsm_output": "X_integrated_totalvi",
    "reference_model_path": "totalvi_model_reference/",
    "query_model_path": "totalvi_model_query/",
}
### VIASH END


def _setup_logger():  # Is there a place where to put common code to not repeat it across components?
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger


def align_proteins_names(adata_reference, mdata_query, adata_query, reference_proteins_key: str, query_proteins_key: str):
    """Make sure that proteins are located in the same .obsm slot in reference and query. Pad query proteins with zeros if they are absent"""
    proteins_reference = adata_reference.obsm[reference_proteins_key]

    # If query has no protein data, put matrix of zeros 
    if not query_proteins_key or query_proteins_key not in mdata_query.mod:
        data = np.zeros((adata_query.n_obs, proteins_reference.shape[1]))
        adata_query.obsm[reference_proteins_key] = DataFrame(columns=proteins_reference.columns, index=adata_query.obs_names, data=data)
    else:
        # Make sure that proteins expression has the same key in query and reference
        adata_query.obsm[reference_proteins_key] = adata_query.obsm[query_proteins_key]

    return adata_query


def convert_mudata_to_anndata(mdata, rna_modality_key, protein_modality_key, input_layer):
    """TOTALVI requires data to be stored in AnnData format with proteins in .obsm slot. This function performs the conversion"""
    adata = mdata.mod[rna_modality_key]

    if protein_modality_key in mdata.mod:
        # Put the proteins modality into .obsm slot
        proteins_reference_adata = mdata.mod[protein_modality_key]

        if input_layer is None:
            proteins = proteins_reference_adata.X
        else:
            proteins = proteins_reference_adata.obsm[input_layer]

        adata.obsm[protein_modality_key] = proteins

    return adata


def build_reference_model(adata_reference) -> scvi.model.TOTALVI:

    vae_reference = scvi.model.TOTALVI(adata_reference, use_layer_norm="both", use_batch_norm="none")
    vae_reference.train()

    return vae_reference

def is_retraining_model() -> bool:
    """Decide, whether reference model should be trained. It happens when no model exists or force_retrain flag is on"""
    
    trained_model_exists = os.path.isdir(par["reference_model_path"]) and ("model.pt" in os.listdir(par["reference_model_path"]))
    return not trained_model_exists or par["force_retrain"]


def map_query_to_reference(mdata_reference, mdata_query, adata_query):
    """Build model on the provided reference if necessary, and map query to the reference"""

    adata_reference = convert_mudata_to_anndata(mdata_reference, rna_modality_key=par["reference_modality"],
                                                protein_modality_key=par["reference_proteins_key"], input_layer=par["input_layer"])

    scvi.model.TOTALVI.setup_anndata(
        adata_reference,
        batch_key=par["obs_batch"],
        protein_expression_obsm_key=par["reference_proteins_key"]
    )

    if is_retraining_model():
        vae_reference = build_reference_model(adata_reference)
    else:
        vae_reference = scvi.model.TOTALVI.load(dir_path=par["reference_model_path"], adata=adata_reference)

    adata_query = align_proteins_names(adata_reference, mdata_query, adata_query, reference_proteins_key=par["reference_proteins_key"],
                                       query_proteins_key=par["query_proteins_key"])

    # Reorder genes and pad missing genes with 0s
    scvi.model.TOTALVI.prepare_query_anndata(adata_query, vae_reference)

    # Train the model for query
    vae_query = scvi.model.TOTALVI.load_query_data(
        adata_query,
        vae_reference
    )
    vae_query.train(200, plan_kwargs=dict(weight_decay=0.0))

    return vae_query, adata_query

def main():
    logger = _setup_logger()

    mdata_query = mudata.read(par["query"].strip())
    adata_query = mdata_query.mod[par["query_modality"]]

    if par["var_input"]:
        logger.info("Subsetting highly variable genes for query")
        selected_genes = adata_query.var_names[adata_query.var[par["var_input"]]]
        adata_query = adata_query[:, selected_genes].copy()

    if par["reference"].endswith(".h5mu"):
        logger.info("Reading reference")
        mdata_reference = mudata.read(par["reference"].strip())

        logger.info("Mapping query to the reference")
        vae_query, adata_query = map_query_to_reference(mdata_reference, mdata_query, adata_query)
    else:
        raise ValueError("Incorrect format of reference, please provide a .h5mu file")

    adata_query.uns["integration_method"] = "totalvi"

    logger.info("Getting the latent representation of query")
    adata_query.obsm[par["obsm_output"]] = vae_query.get_latent_representation()

    mdata_query.mod[par["query_modality"]] = adata_query
    try:
        logger.info("Updating mdata")
        mdata_query.update()
    except KeyError:
        logger.error("Key error was thrown during mdata update. Be careful")

    logger.info("Saving updated query data")
    mdata_query.write_h5mu(par['output'].strip())
    
    logger.info("Saving query model")
    vae_query.save(par["query_model_path"], overwrite=True)

if __name__ == "__main__":
    main()
