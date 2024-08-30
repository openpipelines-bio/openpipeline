import sys
import mudata as mu
import scvi
import numpy as np

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5ad",
    "scvi_reference_model": "resources_test/annotation_test_data/scvi_model.pt",
    "reference_obs_label": "cell_ontology_class",
}
meta = {}
## VIASH END

sys.path.append(meta["resources_dir"])
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

logger.info("Reading the input and reference data")

input_data = mu.read_h5mu(par["input"])
query = input_data.mod[par["modality"]]
reference_data = mu.read_h5mu(par["reference"])
reference = reference_data.mod[par["modality"]]

logger.info(f"Loading the pretrained scVI model from {par['scvi_reference_model']}")
scvi_reference_model = scvi.model.SCVI.load(par["scvi_reference_model"], reference)

logger.info("Setting up scANVI model")

scanvi_ref = scvi.model.SCANVI.from_scvi_model(
    scvi_reference_model,
    unlabeled_category=par["unknown_celltype"],
    labels_key=par["reference_obs_label"],
    )

reference_plan_kwargs = {"lr": par["reference_learning_rate"],
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
    accelerator="auto",
)

logger.info("Updating and training scANVI model with query data")
scvi.model.SCANVI.prepare_query_anndata(query, scanvi_ref, inplace=True)
scanvi_query = scvi.model.SCANVI.load_query_data(query, scanvi_ref)

query_plan_kwargs = {"lr": par["query_learning_rate"],
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
    accelerator="auto",
)

logger.info("Adding latent representation to query data")
query.obsm[par["output_obsm_scanvi_embedding"]] = scanvi_query.get_latent_representation()

logger.info("Running predictions on query data")
query.obs[par["output_obs_predictions"]] = scanvi_query.predict(query)
query.obs[par["output_obs_probability"]] = np.max(scanvi_query.predict(query, soft=True), axis=1)

logger.info("Saving output and model")
input_data.mod[par["modality"]] = query
input_data.write_h5mu(par["output"], compression=par["output_compression"])

if par["output_model"]:
    scanvi_query.save(par["output_model"], overwrite=True)