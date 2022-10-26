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

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "query_modality": "rna",
    "query_proteins_key": None,
    "reference": "HLCA",
    "reference_modality": "rna",
    "reference_proteins_key": None,
    "base_model": "scvi",  # One of ["scvi", "scanvi", "totalvi"]
    "output": "foo.h5mu",
    # scANVI parameters
    "unlabeled_category": "Unknown",
    "labels_key": "leiden",
    "predicted_labels_key": "predicted_labels",
    # Other
    "obsm_output": "X_integrated_{base_model}",  # will be updated with base model
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

PLAN_KWARGS = {
    "reduce_lr_on_plateau": par['reduce_lr_on_plateau'],
    "lr_patience": par['lr_patience'],
    "lr_factor": par['lr_factor'],
}

def build_referrence_model(mdata_reference):
    arches_params = dict(
        use_layer_norm="both",
        use_batch_norm="none")

    adata_reference = mdata_reference.mod[par["reference_modality"]]

    if par["base_model"] == "scvi" or par["base_model"] == "scanvi":
        # Set up the data
        scvi.model.SCVI.setup_anndata(
            adata_reference,
            batch_key=par['obs_batch'],
            layer=par['input_layer']
        )

        # Train the model on the reference
        vae_reference = scvi.model.SCVI(
            adata_reference,
            **arches_params,
            n_hidden=128, #this is the default
            n_latent=30,
            n_layers=2,
            dropout_rate=0.1, #this is the default
            dispersion='gene', #this is the default
            gene_likelihood='nb',
            encode_covariates=True, #Parameterization for better scArches performance -> maybe don't use this always?
            deeply_inject_covariates=False, #Parameterization for better scArches performance -> maybe don't use this always?
            use_observed_lib_size=False, #When size_factors are not passed
        )

        vae_reference.train(
            max_epochs=par['max_epochs'],
            early_stopping=par['early_stopping'],
            early_stopping_monitor=par['early_stopping_monitor'],
            early_stopping_patience=par['early_stopping_patience'],
            early_stopping_min_delta=par['early_stopping_min_delta'],
            plan_kwargs=PLAN_KWARGS,
            check_val_every_n_epoch=1,
            use_gpu=(cuda_is_available() or mps_is_available()),
        )

        if par["base_model"] == "scvi":
            model = scvi.model.SCVI

        elif par["base_model"] == "scanvi":
            vae_reference = scvi.model.SCANVI.from_scvi_model(
                vae_reference,
                unlabeled_category=par["unlabeled_category"],
                labels_key=par["labels_key"],
            )

            vae_reference.train(max_epochs=20, n_samples_per_label=100)
            model = scvi.model.SCANVI

    elif par["base_model"] == "totalvi":

        if par["reference_proteins_key"] in mdata_reference.mod:
            # Put the proteins modality into reference adata
            proteins_reference = mdata_reference.mod[par["reference_proteins_key"]]
            adata_reference.obsm[par["reference_proteins_key"]] = proteins_reference

        scvi.model.TOTALVI.setup_anndata(
            adata_reference,
            batch_key=par['obs_batch'],
            protein_expression_obsm_key=par["reference_proteins_key"]
        )

        vae_reference = scvi.model.TOTALVI(
            adata_reference,
            **arches_params
        )
        vae_reference.train()

        model = scvi.model.TOTALVI

    return adata_reference, vae_reference, model


def align_proteins_names(adata_reference, mdata_query, adata_query):
    proteins_reference = adata_reference.obsm[par["reference_proteins_key"]]

    # If query has no protein data, put matrix of zeros 
    if not par["query_proteins_key"] or par["query_proteins_key"] not in mdata_query.mod:
        data = np.zeros((adata_query.n_obs, proteins_reference.shape[1]))
        adata_query.obsm[par["reference_proteins_key"]] = DataFrame(columns=proteins_reference.columns, index=adata_query.obs_names, data=data)
    else:
        # Make sure that proteins expression has the same key in query and reference
        adata_query.obsm[par["reference_proteins_key"]] = adata_query.obsm[par["query_proteins_key"]]

    return adata_query


def map_query_to_new_reference(mdata_reference, mdata_query, adata_query):
    """Build model on the provided reference and map query to the reference"""

    if par["var_input"]:
        # Subset to HVG
        adata_query = adata_query[:,adata_query.var["var_input"]].copy()

    adata_reference, vae_reference, model = build_referrence_model(mdata_reference)

    if par["base_model"] == "totalvi":
        adata_query = align_proteins_names(adata_reference, mdata_query, adata_query)
        
    # Reorder genes and pad missing genes with 0s
    model.prepare_query_anndata(adata_query, vae_reference)

    # Train the model for query
    vae_query = model.load_query_data(
        adata_query,
        vae_reference
    )
    vae_query.train(200, plan_kwargs=dict(weight_decay=0.0))

    return vae_query, adata_query


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

    print(os.listdir(directory))
    os.system(f"curl  {HLCA_PATH} --output {directory}/HLCA_reference_model.zip")

    print(os.listdir(directory))
    os.system(f"unzip {directory}/HLCA_reference_model.zip -d {directory}")
    print(os.listdir(f"{directory}/HLCA_reference_model"))

    return f"{directory}/HLCA_reference_model"


def main():
    SUPPORTED_BASE_MODELS = set(["scvi", "scanvi", "totalvi"])
    SUPPORTED_REFERENCES = set(["HLCA"])

    if par["reference"] is not None and par["reference"] not in SUPPORTED_REFERENCES:
        raise ValueError(f"{par['reference']} is not supported reference. Please select one of {', '.join(SUPPORTED_REFERENCES)}")

    if par["base_model"] is not None and par["base_model"] not in SUPPORTED_BASE_MODELS:
        raise ValueError(f"{par['base_model']} is not supported base model. Please select one of {', '.join(SUPPORTED_BASE_MODELS)}")

    mdata_query = mudata.read(par["input"].strip())
    adata_query = mdata_query.mod[par["query_modality"]]

    if par["reference"] == "HLCA":
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            model_path = _download_HLCA_reference_model(directory)
            vae_query, adata_query = map_to_existing_reference(adata_query, model_path=model_path)

        par["base_model"] = "scanvi"

    elif par["reference"].endswith(".h5mu"):
        mdata_reference = mudata.read(par["reference"].strip())
        vae_query, adata_query = map_query_to_new_reference(mdata_reference, mdata_query, adata_query)
        
    else:
        raise ValueError(f"Reference {par['reference']} is not supported")

    # Save info about the used model
    adata_query.uns["integration_method"] = par["base_model"]
    output_key = par["obsm_output"].format(base_model=par["base_model"])

    # Get the latent representation
    adata_query.obsm[output_key] = vae_query.get_latent_representation()

    if par["base_model"] == "scanvi":
        adata_query.obs[par["predicted_labels_key"]] = vae_query.predict()

    mdata_query.mod[par["query_modality"]] = adata_query
    mdata_query.write_h5mu(par["output"].strip())

if __name__ == "__main__":
    main()
