import sys
import mudata as mu
import scvi
import numpy as np

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "var_input_gene_names": None,
    "scvi_reference_model": None,
    "scanvi_reference_model": "resources_test/annotation_test_data/scanvi_model",
    "unknown_celltype": "Unkown",
    "output": "output.h5mu",
    "output_obsm_scanvi_embedding": "scanvi_embedding",
    "output_obs_predictions": "scanvi_pred",
    "output_obs_probability": "scanvi_probability",
    "output_model": None,
    "output_compression": None,
    "reference_learning_rate": 1e-3,
    "reference_reduce_lr_on_plateau": True,
    "reference_lr_patience": 25,
    "reference_lr_factor": 0.5,
    "reference_train_size": 0.9,
    "reference_max_epochs": 10,
    "reference_early_stopping": True,
    "reference_early_stopping_patience": 50,
    "query_train_size": 0.9,
    "query_max_epochs": 10,
    "query_learning_rate": 1e-3,
    "query_reduce_lr_on_plateau": True,
    "query_lr_patience": 25,
    "query_lr_factor": 0.5,
    "query_early_stopping": True,
    "query_early_stopping_patience": 50
}
meta = {"resources_dir": "src/annotate/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from query_reference_allignment import set_var_index, cross_check_genes

# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

if (not par["scvi_reference_model"]) and not (par["scanvi_reference_model"]) or (par["scvi_reference_model"] and par["scanvi_reference_model"]):
    raise ValueError("Make sure to provide either an '--scvi_reference_model' or a '--scanvi_reference_model', but not both.")


def main():
    logger.info("Reading the query data")
    # Read in data
    input_data = mu.read_h5mu(par["input"])
    input_modality = input_data.mod[par["modality"]].copy()
    # scANVI requires query and reference gene names to be equivalent 
    input_modality = set_var_index(input_modality, par["var_input_gene_names"])

    if par["scanvi_reference_model"]:

        logger.info(f"Loading the pretrained scANVI model from {par['scanvi_reference_model']} and updating it with the query data {par['input']}")
        scanvi_query = scvi.model.SCANVI.load_query_data(
            input_modality,
            par["scanvi_reference_model"],
            freeze_classifier=True,
            inplace_subset_query_vars=True
            )

    elif par["scvi_reference_model"]:

        logger.info("Reading in the reference model and associated reference data")
        scvi_reference_model = scvi.model.SCVI.load(par["scvi_reference_model"])
        reference = scvi_reference_model.adata


        logger.info("Alligning genes in reference and query dataset")
        # scANVI requires query and reference gene names to be equivalent 
        reference = set_var_index(reference)
        # Subset query dataset based on genes present in reference
        common_ens_ids = cross_check_genes(input_modality.var.index, reference.var.index, min_gene_overlap=par["reference_query_gene_overlap"])
        input_modality = input_modality[:, common_ens_ids]

        logger.info("Instantiating scANVI model from the scVI model")
        scanvi_ref = scvi.model.SCANVI.from_scvi_model(
            scvi_reference_model,
            unlabeled_category=par["unknown_celltype"],
            labels_key=scvi_reference_model.adata_manager._registry["setup_args"]["labels_key"],
            )

        reference_plan_kwargs = {
            "lr": par["reference_learning_rate"],
            "reduce_lr_on_plateau": par['reference_reduce_lr_on_plateau'],
            "lr_patience": par['reference_lr_patience'],
            "lr_factor": par['reference_lr_factor']
            }

        logger.info("Training scANVI model on reference data with celltype labels")

        scanvi_ref.train(
            train_size=par["reference_train_size"],
            max_epochs=par['reference_max_epochs'],
            early_stopping=par['reference_early_stopping'],
            early_stopping_patience=par['reference_early_stopping_patience'],
            plan_kwargs=reference_plan_kwargs,
            check_val_every_n_epoch=1,
            accelerator="auto"
        )

        logger.info(f"Updating scANVI model with query data {par['input']}")
        scvi.model.SCANVI.prepare_query_anndata(input_modality, scanvi_ref, inplace=True)
        scanvi_query = scvi.model.SCANVI.load_query_data(input_modality, scanvi_ref)

    logger.info("Training scANVI model with query data")
    query_plan_kwargs = {
        "lr": par["query_learning_rate"],
        "reduce_lr_on_plateau": par['query_reduce_lr_on_plateau'],
        "lr_patience": par['query_lr_patience'],
        "lr_factor": par['query_lr_factor']
        }

    scanvi_query.train(
        train_size=par["query_train_size"],
        max_epochs=par['query_max_epochs'],
        early_stopping=par['query_early_stopping'],
        early_stopping_patience=par['query_early_stopping_patience'],
        plan_kwargs=query_plan_kwargs,
        check_val_every_n_epoch=1,
        accelerator="auto"
    )

    logger.info("Adding latent representation to query data")
    input_modality.obsm[par["output_obsm_scanvi_embedding"]] = scanvi_query.get_latent_representation()

    logger.info("Running predictions on query data")
    input_modality.obs[par["output_obs_predictions"]] = scanvi_query.predict(input_modality)
    input_modality.obs[par["output_obs_probability"]] = np.max(scanvi_query.predict(input_modality, soft=True), axis=1)

    logger.info("Saving output and model")
    input_data.mod[par["modality"]] = input_modality
    input_data.write_h5mu(par["output"], compression=par["output_compression"])

    if par["output_model"]:
        scanvi_query.save(par["output_model"], overwrite=True)


if __name__ == '__main__':
    main()
