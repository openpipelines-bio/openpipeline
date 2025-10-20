import sys
import mudata as mu
import anndata as ad
import scvi
from scipy.sparse import issparse
from scanpy._utils import check_nonnegative_integers


### VIASH START
par = {
    # inputs
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality_rna": "rna",
    "modality_protein": "prot",
    "input_layer_rna": None,
    "input_layer_protein": None,
    "obs_batch": "sample_id",
    "obs_size_factor": None,
    "obs_categorical_covariate": None,
    "obs_continuous_covariate": None,
    "var_input": None,
    "var_gene_names": None,
    "var_protein_names": None,
    # outputs
    "output": "foo.h5mu",
    "output_compression": "gzip",
    "output_model": "totalvi_model_test/",
    "obsm_output": "X_integrated_totalvi",
    "obsm_normalized_rna_output": "X_totalvi_normalized_rna",
    "obsm_normalized_protein_output": "X_totalvi_normalized_protein",
    # dataset requirements
    "n_obs_min_count": 10,
    "n_var_gene_min_count": 10,
    "n_var_protein_min_count": 3,
    # model arguments
    "n_dimensions_latent_space": 20,
    "gene_dispersion": "gene",
    "protein_dispersion": "protein",
    "gene_likelihood": "nb",
    "latent_distribution": "normal",
    "empirical_protein_background_prior": None,
    "override_missing_proteins": False,
    # training arguments
    "max_epochs": 5,
    "early_stopping": True
}

meta = {"resources_dir": "src/utils/"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from subset_vars import subset_vars
from set_var_index import set_var_index

logger = setup_logger()


def check_validity_anndata(adata, layer, protein_modality, obs_batch, n_obs_min_count, n_var_gene_min_count, n_var_protein_min_count):
    assert check_nonnegative_integers(adata.layers[layer] if layer else adata.X), (
        "Make sure input rna layer contains raw_counts"
    )

    assert len(set(adata.var_names)) == len(adata.var_names), (
        "Dataset contains multiple genes with same gene name."
    )

    # Ensure every obs_batch category has sufficient observations
    assert min(adata.obs[[obs_batch]].value_counts()) > n_obs_min_count, (
        f"Anndata has fewer than {n_obs_min_count} cells."
    )

    assert adata.n_vars > n_var_gene_min_count, (
        f"Anndata has fewer than {n_var_gene_min_count} genes."
    )
    
    assert len(set(adata.var_names)) == len(adata.var_names), (
        "Dataset contains multiple genes with same gene name."
    )

    assert check_nonnegative_integers(adata.obsm[protein_modality]), (
        "Make sure input protein layer contains raw_counts"
    )

    assert adata.obsm[protein_modality].shape[1] > n_var_protein_min_count, (
        f"Anndata has fewer than {n_var_protein_min_count} proteins."
    )


    assert len(set(adata.uns[protein_modality])) == len(adata.uns[protein_modality]), (
        "Dataset contains multiple proteins with same gene name."
    )


def consolidate_modalities_to_anndata(
    mdata: mu.MuData,
    rna_modality: str,
    protein_modality: str,
    input_layer_protein: str | None,
    hvg_var_key: str | None = None,
) -> ad.AnnData:
    """TOTALVI requires data to be stored in AnnData format with protein counts in .obsm slot. This function performs the conversion"""

    adata = mdata.mod[rna_modality].copy()
    if hvg_var_key:
        adata = subset_vars(adata, subset_col=hvg_var_key)

    adata = set_var_index(adata, var_name=par["var_gene_names"])

    # Put the proteins modality into .obsm slot
    prot_adata = mdata.mod[protein_modality].copy()
    set_var_index(prot_adata, var_name=par["var_protein_names"])

    protein_layer = prot_adata.layers[input_layer_protein] if input_layer_protein else prot_adata.X

    if issparse(protein_layer):
        protein_layer = protein_layer.toarray()

    adata.obsm[protein_modality] = protein_layer
    adata.uns[protein_modality] = prot_adata.var_names

    return adata


def main():
    mdata = mu.read_h5mu(par["input"])
    
    logger.info("Preparing data for TOTALVI...")
    adata = consolidate_modalities_to_anndata(
        mdata,
        par["modality_rna"],
        par["modality_protein"],
        par["input_layer_protein"],
        par["var_input"],
    )

    check_validity_anndata(
        adata,
        par["input_layer_rna"],
        par["modality_protein"],
        par["obs_batch"],
        par["n_obs_min_count"],
        par["n_var_gene_min_count"],
        par["n_var_protein_min_count"]
    )

    logger.info("Loading data for TOTALVI...")
    scvi.model.TOTALVI.setup_anndata(
        adata,
        batch_key=par["obs_batch"],
        layer=par["input_layer_rna"],
        protein_expression_obsm_key=par["modality_protein"],
        protein_names_uns_key=par["modality_protein"],
        size_factor_key=par["obs_size_factor"],
        categorical_covariate_keys=par["obs_categorical_covariate"],
        continuous_covariate_keys=par["obs_continuous_covariate"],
    )

    logger.info("Preparing TOTALVI model...")
    model = scvi.model.TOTALVI(
        adata,
        n_latent=par["n_dimensions_latent_space"],
        gene_dispersion=par["gene_dispersion"],
        protein_dispersion=par["protein_dispersion"],
        gene_likelihood=par["gene_likelihood"],
        latent_distribution=par["latent_distribution"],
        empirical_protein_background_prior=par["empirical_protein_background_prior"],
        override_missing_proteins=par["override_missing_proteins"],
    )
    
    logger.info("Training TOTALVI model...")
    model.train(
        max_epochs=par["max_epochs"],
        early_stopping=par["early_stopping"],
        check_val_every_n_epoch=1,
        accelerator="auto",
    )
    
    logger.info("Getting the latent representation...")
    mdata.mod[par["modality_rna"]].obsm[par["obsm_output"]] = model.get_latent_representation()
    
    logger.info("Getting normalized expression...")
    norm_rna, norm_protein = model.get_normalized_expression()
    mdata.mod[par["modality_rna"]].obsm[par["obsm_normalized_rna_output"]] = norm_rna.to_numpy()
    if par["modality_protein"] in mdata.mod:
        mdata.mod[par["modality_protein"]].obsm[par["obsm_normalized_protein_output"]] = norm_protein.to_numpy()

    logger.info("Saving integrated data...")
    mdata.write_h5mu(par["output"], compression=par["output_compression"])
    
    if par["output_model"]:
        logger.info("Saving TOTALVI model...")
        model.save(par["output_model"], overwrite=True)


if __name__ == "__main__":
    main()
