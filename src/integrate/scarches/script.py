import logging
from pyexpat import model
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
    # Other
    "obsm_output": "X_integrated_{model_name}",
    "obs_batch": "sample_id",
    "input_layer": None,
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

def _read_model_name_from_registry(model_path) -> str:
    """Read registry with information about the model, return the model name"""
    registry = scvi.model.base.BaseModelClass.load_registry(model_path)
    return registry["model_name"]


def _detect_base_model(model_path) -> scvi.model.base.BaseModelClass:
    """Read from the model's file which scvi_tools model it contains"""

    names_to_models_map = {
        "AUTOZI": scvi.model.AUTOZI,
        "CondSCVI": scvi.model.CondSCVI,
        "DestVI": scvi.model.DestVI,
        "LinearSCVI": scvi.model.LinearSCVI,
        "PEAKVI": scvi.model.PEAKVI,
        "SCANVI": scvi.model.SCANVI,
        "SCVI": scvi.model.SCVI,
        "TOTALVI": scvi.model.TOTALVI,
        "MULTIVI": scvi.model.MULTIVI,
        "AmortizedLDA": scvi.model.AmortizedLDA,
        "JaxSCVI": scvi.model.JaxSCVI,
    }
    
    return names_to_models_map[_read_model_name_from_registry(model_path)]


def map_to_existing_reference(adata_query, model_path, check_val_every_n_epoch=1):
    """
    A function to map the query data to the reference atlas

    Input:
        * adata_query: An AnnData object with the query
        * model_path: The reference model directory
    
    Output:
        * vae_query: the trained scvi_tools model
        * adata_query: The AnnData object with the query preprocessed for the mapping to the reference
    """
    model = _detect_base_model(model_path)
    model.prepare_query_anndata(adata_query, model_path)

    # Load query data into the model
    vae_query = model.load_query_data(
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
            use_gpu=(cuda_is_available() or mps_is_available())
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
            model_name = _read_model_name_from_registry(directory)
        
    else:
        raise ValueError(f"Reference {par['reference']} is not supported")

    # Save info about the used model
    adata_query.uns["integration_method"] = model_name

    logger.info("Trying to write latent representation")
    output_key = par["obsm_output"].format(model_name=model_name)
    adata_query.obsm[output_key] = vae_query.get_latent_representation()

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
