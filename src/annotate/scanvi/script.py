import sys
import mudata as mu
import anndata as ad
import scvi

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

logger.info("Reading input data and SCVI model")

input_data = mu.read_h5mu(par["input"])
query = input_data.mod[par["modality"]]
reference = ad.read_h5ad(par["reference"])
scvi_reference_model = scvi.model.SCVI.load(par["scvi_reference_model"], reference)

logger.info("Setting up scANVI model")

scanvi_ref = scvi.model.SCANVI.from_scvi_model(
    scvi_reference_model,
    unlabeled_category="Unknown",
    labels_key=par["reference_obs_label"],
    )

plan_kwargs = {"lr": par["learning_rate"],
               "reduce_lr_on_plateau": par['reduce_lr_on_plateau'],
               "lr_patience": par['lr_patience'],
               "lr_factor": par['lr_factor']
               }

scanvi_ref.train(
    train_size=par["train_size"],
    max_epochs=par['max_epochs'],
    early_stopping=par['early_stopping'],
    early_stopping_patience=par['early_stopping_patience'],
    plan_kwargs=plan_kwargs,
    check_val_every_n_epoch=1,
    accelerator="auto",
)

SCANVI_LATENT_KEY = "X_scANVI"
reference.obsm[SCANVI_LATENT_KEY] = scanvi_ref.get_latent_representation()

logger.info("Updating with query data")

scvi.model.SCANVI.prepare_query_anndata(query, scanvi_ref, inplace=True)
scanvi_query = scvi.model.SCANVI.load_query_data(query, scanvi_ref)


scanvi_query.train(
    train_size=par["train_size"],
    max_epochs=par['max_epochs'],
    early_stopping=par['early_stopping'],
    early_stopping_patience=par['early_stopping_patience'],
    plan_kwargs=plan_kwargs,
    check_val_every_n_epoch=1,
    accelerator="auto",
)

logger.info("Running prediction")

query.obsm[SCANVI_LATENT_KEY] = scanvi_query.get_latent_representation()
query.obs[par["input_obs_label"]] = scanvi_query.predict()

logger.info("Saving output and model")

input_data.mod[par["modality"]] = query
input_data.write_h5mu(par["output"], compression=par["output_compression"])

if par["output_model"]:
    scanvi_query.save(par["output_model"], overwrite=True)