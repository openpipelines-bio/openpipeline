import logging
import mudata
import scvi
from torch.cuda import is_available as cuda_is_available
try:
    from torch.backends.mps import is_available as mps_is_available
except ModuleNotFoundError:
    # Older pytorch versions
    # MacOS GPUs
    def mps_is_available():
        return False

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "output": "foo.h5mu",
    # scANVI parameters
    "unlabeled_category": "Unknown",
    "labels_key": "leiden",
    "predicted_labels_key": "predicted_labels",
    # Other
    "obsm_output": "X_integrated_scanvi",
    "var_input": None,
    "obs_batch": "sample_id",
    "input_layer": None,
    "reduce_lr_on_plateau": True,
    "lr_factor": 0.6,
    "lr_patience": 30,
    "early_stopping": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 0,
    "max_epochs": 500}
### VIASH END

def _setup_logger():
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger

PLAN_KWARGS = {
    "reduce_lr_on_plateau": par['reduce_lr_on_plateau'],
    "lr_patience": par['lr_patience'],
    "lr_factor": par['lr_factor'],
}


def map_to_existing_reference(adata_query, model_path, check_val_every_n_epoch=1) -> np.ndarray:
    """
    A function to map the query data to the reference atlas

    Input:
        * adata_query: An AnnData object with the query
        * model_path: The reference model directory

    Output:
        A numpy.ndarray containing the embedding coordinates of the query data in the HLCA latent space.
        If `output_model` is set to True, then the trained reference model is also output

    """
    print("SCVI-tools version:", scvi.__version__)

    scvi.model.SCANVI.prepare_query_anndata(adata_query, model_path)

    # Load query data into the model
    vae_query = scvi.model.SCANVI.load_query_data(
            adata_query,
            model_path,
            freeze_dropout=True
    )

    # Train scArches model for query mapping
    vae_query.train(
            max_epochs=par["max_epochs"],
            early_stopping=par['early_stopping'],
            early_stopping_monitor=par['early_stopping_monitor'],
            early_stopping_patience=par['early_stopping_patience'],
            early_stopping_min_delta=par['early_stopping_min_delta'],
            check_val_every_n_epoch=check_val_every_n_epoch,
            use_gpu=(cuda_is_available() or mps_is_available()),
            plan_kwargs=PLAN_KWARGS
    )

    return vae_query, adata_query


def _download_HLCA_reference_model(directory):
    import os

    HLCA_PATH = "https://zenodo.org/record/6337966/files/HLCA_reference_model.zip"

    os.system(f"curl  {HLCA_PATH} --output {directory}/HLCA_reference_model.zip")
    os.system(f"unzip {directory}/HLCA_reference_model.zip -d {directory}")

    scvi.model.SCANVI.convert_legacy_save(f"{directory}/HLCA_reference_model", f"{directory}/HLCA_reference_model_new")

    return f"{directory}/HLCA_reference_model_new"


def main():
    SUPPORTED_REFERENCES = ["HLCA"]

    if par["reference"] is not None and par["reference"] not in SUPPORTED_REFERENCES:
        raise ValueError(f"{par['reference']} is not supported reference. Please select one of {', '.join(SUPPORTED_REFERENCES)}")

    logger = _setup_logger()

    mdata_query = mudata.read(par["input"].strip())
    adata_query = mdata_query.mod[par["modality"]]

    if par["reference"] == "HLCA":
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            model_path = _download_HLCA_reference_model(directory)
            vae_query, adata_query = map_to_existing_reference(adata_query, model_path=model_path)
        
    else:
        raise ValueError(f"Reference {par['reference']} is not supported")

    # Save info about the used model
    adata_query.uns["integration_method"] = "scanvi"

    logger.info("Trying to write latent representation")
    adata_query.obsm[par["obsm_output"]] = vae_query.get_latent_representation()

    logger.info("Predicting labels for SCANVI")
    adata_query.obs[par["predicted_labels_key"]] = vae_query.predict()

    logger.info("Reassigning modality")
    mdata_query.mod[par["modality"]] = adata_query
    try:
        mdata_query.update()
    except KeyError:
        logger.error("Key error was thrown during mdata update. Be careful")

    logger.info("Saving h5mu file")
    mdata_query.write_h5mu(par["output"].strip())

if __name__ == "__main__":
    main()
