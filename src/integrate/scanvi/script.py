
import mudata
import scvi

### VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "scvi_model": "resources_test/annotation_test_data/scvi_model",
    "obs_labels": "cell_ontology_class",
    "unlabeled_category": "Unknown",
    "modality": "rna",
    "input_layer": None,
    "var_input": None,
    "output": "foo.h5mu",
    "var_gene_names": "ensemblid",
    "obsm_output": "X_scvi_integrated",
    "early_stopping": True,
    "early_stopping_monitor": "elbo_validation",
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 0,
    "reduce_lr_on_plateau": True,
    "lr_factor": 0.6,
    "lr_patience": 30,
    "max_epochs": 5,
    "n_obs_min_count": 10,
    "n_var_min_count": 10,
    "output_model": "test/",
    "output_compression": "gzip",
    "output_model": "resources_test/annotation_test_data/scanvi_model"
}

meta = {"resources_dir": "src/utils"}
### VIASH END

import sys

sys.path.append(meta["resources_dir"])

from subset_vars import subset_vars
from compress_h5mu import write_h5ad_to_h5mu_with_compression
from setup_logger import setup_logger
from set_var_index import set_var_index

logger = setup_logger()

def main():
    logger.info("Reading input data...")
    adata = mudata.read_h5ad(par["input"].strip(), mod=par["modality"])

    if par["var_input"]:
        # Subset to HVG
        adata_subset = subset_vars(adata, subset_col=par["var_input"]).copy()
    else:
        adata_subset = adata.copy()

    # Sanitize gene names and set as index of the AnnData object
    adata_subset = set_var_index(adata_subset, par["var_gene_names"])
    
    logger.info(f"Loading pre-trained scVI model from {par['scvi_model']}")
    scvi_model = scvi.model.SCVI.load(
        par["scvi_model"],
        adata_subset,
        accelerator="auto",
        device="auto",
    )
    
    logger.info("Instantiating scANVI model from scVI model...")
    vae_uns = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        unlabeled_category=par["unlabeled_category"],
        labels_key=par["obs_labels"],
        adata=adata_subset
    )

    plan_kwargs = {
        "reduce_lr_on_plateau": par["reduce_lr_on_plateau"],
        "lr_patience": par["lr_patience"],
        "lr_factor": par["lr_factor"],
    }

    logger.info("Training scANVI model...")
    # Train the model
    vae_uns.train(
        max_epochs=par["max_epochs"],
        early_stopping=par["early_stopping"],
        early_stopping_monitor=par["early_stopping_monitor"],
        early_stopping_patience=par["early_stopping_patience"],
        early_stopping_min_delta=par["early_stopping_min_delta"],
        plan_kwargs=plan_kwargs,
        check_val_every_n_epoch=1,
        accelerator="auto",
    )
    # Note: train_size=1.0 should give better results, but then can't do early_stopping on validation set

    logger.info("Performing scANVI integration...")
    # Get the latent output
    adata.obsm[par["obsm_output"]] = vae_uns.get_latent_representation()

    logger.info("Writing output data...")
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], adata, par["output_compression"]
    )
    if par["output_model"]:
        vae_uns.save(par["output_model"], overwrite=True)


if __name__ == "__main__":
    main()
