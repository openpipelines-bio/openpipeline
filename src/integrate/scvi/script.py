from scanpy._utils import check_nonnegative_integers
import mudata
import scvi
import sys

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_layer": None,
    "obs_batch": "sample_id",
    "obs_labels": None,
    "var_input": None,
    "output": "foo.h5mu",
    "obsm_output": "X_scvi_integrated",
    "obs_categorical_covariate": None,
    "obs_continuous_covariate": None,
    "n_hidden_nodes": 128,
    "n_dimensions_latent_space": 30,
    "n_hidden_layers": 2,
    "dropout_rate": 0.1,
    "dispersion": "gene",
    "gene_likelihood": "nb",
    "use_layer_normalization": "both",
    "use_batch_normalization": "none",
    "encode_covariates": True,
    "deeply_inject_covariates": False,
    "use_observed_lib_size": False,
    "early_stopping": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 0,
    "reduce_lr_on_plateau": True,
    "lr_factor": 0.6,
    "lr_patience": 30,
    "max_epochs": 50,
    "n_obs_min_count": 10,
    "n_var_min_count": 10,
    "output_model": "scvi/",
    "obs_size_factor": None,
    "output_compression": "gzip",
    "var_input_gene_names": None,
    "scvi_reference_model": "scvi",
    "unknown_celltype": "Unknown",
    "input_reference_gene_overlap": 100,
}

meta = {"resources_dir": "src/utils"}
## VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from subset_vars import subset_vars
from set_var_index import set_var_index

logger = setup_logger()


# TODO: optionally, move to qa
# https://github.com/openpipelines-bio/openpipeline/issues/435
def check_validity_anndata(adata, layer, obs_batch, n_obs_min_count, n_var_min_count):
    assert check_nonnegative_integers(
        adata.layers[layer] if layer else adata.X
    ), "Make sure input adata contains raw_counts"

    assert len(set(adata.var_names)) == len(
        adata.var_names
    ), "Dataset contains multiple genes with same gene name."

    # Ensure every obs_batch category has sufficient observations
    assert (
        min(adata.obs[[obs_batch]].value_counts()) > n_obs_min_count
    ), f"Anndata has fewer than {n_obs_min_count} cells."

    assert (
        adata.n_vars > n_var_min_count
    ), f"Anndata has fewer than {n_var_min_count} genes."


def align_query_with_reference(query, par, reference_registry):
    setup_args = reference_registry["setup_args"]
    # Allign layers
    if not setup_args["layer"] == par["input_layer"]:
        logger.info(
            f"Reference model layer {setup_args['layer']} is different from query layer {par['input_layer']}."
        )
        logger.info("Aligning layers.")
        # Can be overwritten since final results will be copied back to original input data
        query.layers[setup_args["layer"]] = query.X

    # Allign batch key
    if not setup_args["batch_key"] == par["obs_batch"]:
        logger.info(
            f"Reference model batch key {setup_args['batch_key']} is different from query batch key {par['obs_batch']}."
        )
        logger.info("Aligning batch keys.")
        # Can be overwritten since final results will be copied back to original input data
        query.obs[setup_args["batch_key"]] = query.obs[par["obs_batch"]]

    if not setup_args["labels_key"] == par["obs_labels"]:
        # Can be overwritten since final results will be copied back to original input data
        if par["obs_labels"]:
            logger.info(
                f"Reference model labels key {setup_args['labels_key']} is different from query labels key {par['obs_labels']}."
            )
            logger.info("Aligning labels keys.")
            query.obs[setup_args["labels_key"]] = query.obs[par["obs_labels"]]
        else:
            logger.info(
                f"No labels key provided for query data. Setting all labels to {par['unknown_celltype']}."
            )
            query.obs[setup_args["labels_key"]] = par["unknown_celltype"]

    return query


def main():
    mdata = mudata.read(par["input"].strip())
    adata = mdata.mod[par["modality"]]
    input_modality = adata.copy()
    # scVI requires query and reference gene names to be equivalent
    input_modality = set_var_index(input_modality, par["var_input_gene_names"])

    if par["scvi_reference_model"]:
        logger.info(
            f"Loading the pretrained scVI model from {par['scvi_reference_model']}"
        )
        scvi_reference_model = scvi.model.SCVI.load(par["scvi_reference_model"])

        logger.info("Alligning genes in reference and query dataset")

        # Mapping of fields
        scvi.model.SCVI.setup_anndata(
            input_modality,
            batch_key=par["obs_batch"],
            layer=par["input_layer"],
            labels_key=par["obs_labels"],
            size_factor_key=par["obs_size_factor"],
            categorical_covariate_keys=par["obs_categorical_covariate"],
            continuous_covariate_keys=par["obs_continuous_covariate"],
        )

        # Allignment of features
        # .varm an be cleared since final results will be copied back to original input data
        # Avoids ValueErrors when features of query differ from reference
        input_modality.varm.clear()
        scvi.model.SCVI.prepare_query_anndata(input_modality, scvi_reference_model)

        input_modality = align_query_with_reference(
            input_modality, par, scvi_reference_model.registry_
        )

        logger.info("Updating pretrained scVI reference model with query data")
        vae_uns = scvi.model.SCVI.load_query_data(input_modality, scvi_reference_model)

    else:
        if par["var_input"]:
            # Subset to HVG
            input_modality = subset_vars(input_modality, subset_col=par["var_input"])

        check_validity_anndata(
            input_modality,
            par["input_layer"],
            par["obs_batch"],
            par["n_obs_min_count"],
            par["n_var_min_count"],
        )
        # Set up the data
        scvi.model.SCVI.setup_anndata(
            input_modality,
            batch_key=par["obs_batch"],
            layer=par["input_layer"],
            labels_key=par["obs_labels"],
            size_factor_key=par["obs_size_factor"],
            categorical_covariate_keys=par["obs_categorical_covariate"],
            continuous_covariate_keys=par["obs_continuous_covariate"],
        )

        # Set up the model
        vae_uns = scvi.model.SCVI(
            input_modality,
            n_hidden=par["n_hidden_nodes"],
            n_latent=par["n_dimensions_latent_space"],
            n_layers=par["n_hidden_layers"],
            dropout_rate=par["dropout_rate"],
            dispersion=par["dispersion"],
            gene_likelihood=par["gene_likelihood"],
            use_layer_norm=par["use_layer_normalization"],
            use_batch_norm=par["use_batch_normalization"],
            encode_covariates=par[
                "encode_covariates"
            ],  # Default (True) is for better scArches performance -> maybe don't use this always?
            deeply_inject_covariates=par[
                "deeply_inject_covariates"
            ],  # Default (False) for better scArches performance -> maybe don't use this always?
            use_observed_lib_size=par[
                "use_observed_lib_size"
            ],  # When size_factors are not passed
        )

    plan_kwargs = {
        "reduce_lr_on_plateau": par["reduce_lr_on_plateau"],
        "lr_patience": par["lr_patience"],
        "lr_factor": par["lr_factor"],
    }

    # Train the model
    vae_uns.train(
        max_epochs=par["max_epochs"],
        early_stopping=par["early_stopping"],
        early_stopping_monitor=par["early_stopping_monitor"],
        early_stopping_patience=par["early_stopping_patience"],
        early_stopping_min_delta=par["early_stopping_min_delta"],
        plan_kwargs=plan_kwargs,
        check_val_every_n_epoch=1,
        accelerator="cpu",
        devices="auto",
    )
    # Note: train_size=1.0 should give better results, but then can't do early_stopping on validation set

    # Get the latent output
    adata.obsm[par["obsm_output"]] = vae_uns.get_latent_representation()

    mdata.mod[par["modality"]] = adata
    mdata.write_h5mu(par["output"].strip(), compression=par["output_compression"])
    if par["output_model"]:
        vae_uns.save(par["output_model"], overwrite=True, save_anndata=True)


if __name__ == "__main__":
    main()
