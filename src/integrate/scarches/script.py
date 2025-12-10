import sys
import mudata
from anndata import AnnData
import scvi
import pandas as pd
import numpy as np
import torch
from tempfile import TemporaryDirectory
from pathlib import Path
import pickle

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "layer": None,
    "input_obs_batch": "sample_id",
    "input_obs_label": None,
    "input_var_gene_names": None,
    "input_obs_categorical_covariate": None,
    "input_obs_continuous_covariate": None,
    "unknown_celltype_label": "Unknown",
    "input_obs_size_factor": None,
    "reference": "resources_test/annotation_test_data/scanvi_model",
    "modality": "rna",
    "output": "foo.h5mu",
    "model_output": "test",
    # Other
    "obsm_output": "X_integrated_scanvi",
    "obs_output_probabilities": "scanvi_proba",
    "obs_output_predictions": "scanvi_pred",
    "early_stopping": None,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 0,
    "max_epochs": 10,
    "output_compression": "gzip",
    "reference_class": "SCANVI",
}
meta = {"resources_dir": "src/utils", "temp_dir": "/tmp"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from set_var_index import set_var_index
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def _read_model_name_from_registry(model_path) -> str:
    """Read registry with information about the model, return the model name"""
    registry = scvi.model.base.BaseModelClass.load_registry(model_path)
    return registry["model_name"]


def _detect_base_model(model_path):
    """Read from the model's file which scvi_tools model it contains"""
    return getattr(scvi.model, _read_model_name_from_registry(model_path))


def _validate_obs_metadata_params(model_registry, model_name):
    """
    Validates .obs metadata par variables that are required by scvi-tools models.

    This function checks that necessary .obs metadata field names (batch, labels, size factors)
    specified in the model registry have the required corresponding parameters provided in the input.

    Parameters
    ----------
    model_registry : dict
        Dictionary containing model configuration and requirements.
    model_name : str
        Name of the model being validated.
    """

    # Maps the keys of the model registry to their corresponding par keys containing the .obs metadata
    registry_par_mapper = {
        "batch_key": "input_obs_batch",
        "labels_key": "input_obs_label",
        "size_factor_key": "input_obs_size_factor",
        "categorical_covariate_keys": "input_obs_categorical_covariate",
        "continuous_covariate_keys": "input_obs_continuous_covariate",
    }

    for registry_key, key in registry_par_mapper.items():
        model_registry_key = model_registry.get(registry_key)
        par_key = par[key]

        if model_registry_key and not par_key:
            if key != "input_obs_label" or not model_registry.get("unlabeled_category"):
                raise ValueError(
                    f"The provided {model_name} model requires `--{key}` to be provided."
                )

        elif par_key and not model_registry_key:
            logger.warning(
                f"`--{key}` was provided but is not used in the provided {model_name} model."
            )

        elif model_registry_key and par_key:
            if key in [
                "input_obs_categorical_covariate",
                "input_obs_continuous_covariate",
            ]:
                if len(par_key) != len(model_registry_key):
                    raise ValueError(
                        f"The number of provided covariates in `--{key}` ({par_key}) does not match the number of covariates used in the provided model ({model_registry_key})."
                    )


def _align_query_with_registry(adata_query, model_path):
    """
    Creates a qeury AnnData object with the expected structure and metadata fields that are aligned with the pre-trained reference model.

    Parameters
    ----------
    adata_query : AnnData
        The query AnnData object to be aligned with the model structure.
    model_path : str
        Path to the directory containing the pre-trained model.

    Returns
    -------
    AnnData
        A new AnnData object with structure and metadata aligned to match the
        requirements of the pre-trained model.
    """

    model = _detect_base_model(model_path)
    model_name = _read_model_name_from_registry(model_path)
    model_registry = model.load_registry(model_path)["setup_args"]
    _validate_obs_metadata_params(model_registry, model_name)

    # Sanitize gene names and set as index of the AnnData object
    # all scArches VAE models expect gene names to be in the .var index
    adata_query = set_var_index(
        adata_query, par["input_var_gene_names"], par["sanitize_ensembl_ids"]
    )

    # align layer
    query_layer = (
        adata_query.X if not par["layer"] else adata_query.layers[par["layer"]]
    )
    var_index = adata_query.var_names.tolist()
    obs_index = adata_query.obs_names.tolist()

    # align observations
    query_obs = {}

    ## batch_key, size_factor_key
    simple_mappings = {
        "batch_key": "input_obs_batch",  # relevant for AUTOZI, LinearSCVI, PEAKVI, SCANVI, SCVI, TOTALVI, MULTIVI, JaxSCVI
        "size_factor_key": "input_obs_size_factor",  # relevant for SCANVI, SCVI, TOTALVI, MULTIVI
    }

    for registry_key, par_key in simple_mappings.items():
        if model_registry.get(registry_key):
            query_obs[model_registry[registry_key]] = adata_query.obs[
                par[par_key]
            ].tolist()

    ## labels-key, relevant for AUTOZI, CondSCVI, LinearSCVI, PEAKVI, SCANVI, SCVI
    if model_registry.get("labels_key"):
        if par["input_obs_label"]:
            query_obs[model_registry["labels_key"]] = adata_query.obs[
                par["input_obs_label"]
            ].tolist()
        else:
            adata_query.obs[model_registry["labels_key"]] = model_registry[
                "unlabeled_category"
            ]
            query_obs[model_registry["labels_key"]] = adata_query.obs[
                model_registry["labels_key"]
            ].tolist()

    ## covariate keys
    covariate_mappings = {
        "categorical_covariate_keys": "input_obs_categorical_covariate",
        "continuous_covariate_keys": "input_obs_continuous_covariate",
    }
    for registry_key, par_key in covariate_mappings.items():
        if model_registry.get(registry_key):
            for covariate, par_covariate in zip(
                model_registry[registry_key], par[par_key]
            ):
                query_obs[covariate] = adata_query.obs[par_covariate].tolist()

    obs = pd.DataFrame(query_obs, index=obs_index)
    var = pd.DataFrame(index=var_index)

    aligned_query_anndata = AnnData(X=query_layer, obs=obs, var=var)
    if model_registry["layer"]:
        aligned_query_anndata.layers[model_registry["layer"]] = query_layer

    return aligned_query_anndata


def map_to_existing_reference(adata_query, model_path, check_val_every_n_epoch=1):
    """
    A function to map the query data to the reference atlas.

    Input:
        * adata_query: An AnnData object with the query
        * model_path: The reference model directory

    Output:
        * vae_query: the trained scvi_tools model
        * adata_query: The AnnData object with the query preprocessed for the mapping to the reference
    """
    model = _detect_base_model(model_path)

    # Keys of the AnnData query object need to match the exact keys in the reference model registry

    aligned_adata_query = _align_query_with_registry(adata_query, model_path)

    try:
        model.prepare_query_anndata(aligned_adata_query, model_path)
    except ValueError:
        logger.warning(
            "ValueError thrown when preparing adata for mapping. Clearing .varm field to prevent it"
        )
        aligned_adata_query.varm.clear()
        try:
            model.prepare_query_anndata(aligned_adata_query, model_path)
        except ValueError:
            raise ValueError(
                f"Could not perform model.prepare_model_anndata, likely because the model was trained with different var names then were found in the index. \n\nmodel var_names: {model.prepare_query_anndata(aligned_adata_query, model_path, return_reference_var_names=True).tolist()} \n\nquery data var_names: {aligned_adata_query.var_names.tolist()}  "
            )

    try:
        # Load query data into the model
        vae_query = model.load_query_data(
            aligned_adata_query, model_path, freeze_dropout=True
        )
    except KeyError:
        # Older models do not have the 'setup_method_name' key saved with the model.
        # Assume these were generated from an anndata object.
        model_torch = torch.load(
            f"{model_path}/model.pt", weights_only=False, map_location="cpu"
        )
        model_torch["attr_dict"]["registry_"]["setup_method_name"] = "setup_anndata"

        with TemporaryDirectory(dir=meta["temp_dir"]) as tempdir:
            temp_file_name = Path(tempdir) / "model.pt"
            torch.save(model_torch, temp_file_name)
            del model_torch
            vae_query = model.load_query_data(
                aligned_adata_query, tempdir, freeze_dropout=True
            )

    # Train scArches model for query mapping
    vae_query.train(
        max_epochs=par["max_epochs"],
        early_stopping=par["early_stopping"],
        early_stopping_monitor=par["early_stopping_monitor"],
        early_stopping_patience=par["early_stopping_patience"],
        early_stopping_min_delta=par["early_stopping_min_delta"],
        check_val_every_n_epoch=check_val_every_n_epoch,
        accelerator="auto",
    )

    return vae_query


def _convert_object_dtypes_to_strings(adata):
    """Convert object dtypes in .var and .obs to string to prevent error when saving file"""

    def convert_cols(df):
        object_cols = df.columns[df.dtypes == "object"]
        for col in object_cols:
            df[col] = df[col].astype(str)
        return df

    adata.var = convert_cols(adata.var)
    adata.obs = convert_cols(adata.obs)

    return adata


def _get_model_path(model_path: str):
    """Obtain path to the directory with reference model. If the proposed `model_path` is a .zip archive, unzip it. If nesessary, convert model to the new format

    Parameters
    ----------
    model_path : str
        Path to a directory, where to search for the model or to a zip file containing the model

    Returns
    -------
    Path to a directory with reference model in format of scvi-tools>=0.15
    """
    import os
    import zipfile
    import tempfile
    from pathlib import Path

    if os.path.isdir(model_path) and "model.pt" in os.listdir(model_path):
        # Probably, the `model_path` already contains model in the output format of scvi-tools>=0.15
        return model_path

    # The model either has old format or is a zip file downloaded from Zenodo
    new_directory = Path(tempfile.TemporaryDirectory().name)

    if zipfile.is_zipfile(model_path):
        with zipfile.ZipFile(model_path) as archive:
            archive.extractall(new_directory)
            model_dir = next(new_directory.glob("**/*.pt")).parent

    else:
        model_dir = next(Path(model_path).glob("**/*.pt")).parent

    if "model_params.pt" in os.listdir(model_dir):
        try:
            # The current approach is to store the class with the model in a registry
            model_class = _detect_base_model(model_dir)
        except ValueError:
            # Failed to load the registry, try
            with (model_dir / "attr.pkl").open("rb") as open_file:
                attrs = pickle.load(open_file)
            try:
                model_name = attrs["registry_"]["model_name"]
            except KeyError:
                # With older scvi-tools models, the model name was not stored with the model...
                model_name = par["reference_class"]
                if not model_name:
                    raise ValueError(
                        "Could not load reference model. This might be because it is a legacy model for which the model type was not saved with the model. Please provide it manually using the 'reference_class' argument."
                    )
                if not hasattr(scvi.model, model_name):
                    raise ValueError(
                        f"Requested to load legacy model using type {model_name}, but such a class does not exist in `scvi.models`."
                    )
            model_class = getattr(scvi.model, model_name)

        # The model is in the `directory`, but it was generated with scvi-tools<0.15
        # TODO: for new references (that could not be SCANVI based), we need to check the base class somehow. Reading registry does not work with models generated by scvi-tools<0.15
        # Here I assume that the reference model is for HLCA and thus is SCANVI based
        converted_model_path = os.path.join(model_dir, "converted")
        model_class.convert_legacy_save(model_dir, converted_model_path)
        return converted_model_path

    elif "model.pt" in os.listdir(model_dir):
        # Archive contained model in the new format, so just return the directory
        return model_dir

    else:
        raise ValueError(
            "Cannot find model in the provided reference path. Please, provide a path or a link to the directory with reference model. For HLCA use https://zenodo.org/record/6337966/files/HLCA_reference_model.zip"
        )


def main():
    logger.info("Reading %s, modality %s", par["input"], par["modality"])
    adata = mudata.read_h5ad(par["input"].strip(), mod=par["modality"])
    adata_query = adata.copy()

    model_path = _get_model_path(par["reference"])
    vae_query = map_to_existing_reference(adata_query, model_path=model_path)
    model_name = _read_model_name_from_registry(model_path)

    # Save info about the used model
    adata.uns["integration_method"] = model_name

    logger.info("Trying to write latent representation")
    output_key = par["obsm_output"].format(model_name=model_name)
    adata.obsm[output_key] = vae_query.get_latent_representation()

    # In addition to integrating the data, scANVI also performs label annotations
    if model_name == "SCANVI":
        adata.obs[par["obs_output_predictions"]] = vae_query.predict()
        adata.obs[par["obs_output_probabilities"]] = np.max(
            vae_query.predict(soft=True), axis=1
        )

    logger.info("Converting dtypes")
    adata = _convert_object_dtypes_to_strings(adata)

    logger.info(
        "Saving h5mu file to %s with compression %s",
        par["output"],
        par["output_compression"],
    )
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], adata, par["output_compression"]
    )

    logger.info("Saving model")
    vae_query.save(par["model_output"], overwrite=True)


if __name__ == "__main__":
    main()
