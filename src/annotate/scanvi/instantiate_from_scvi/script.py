import sys
import scvi

## VIASH START
par = {
    "input": "resources_test/annotation_test_data/scvi_model",
    "output": "resources_test/annotation_test_data/scanvi_model",
    "learning_rate": 1e-3,
    "reduce_lr_on_plateau": True,
    "lr_patience": 25,
    "lr_factor": 0.5,
    "train_size": 0.9,
    "max_epochs": 10,
    "early_stopping": True,
    "early_stopping_patience": 50,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

def main():

    logger.info("Reading in the SCVI model")
    scvi_reference_model = scvi.model.SCVI.load(par["input"])

    logger.info("Instantiating scANVI model from the scVI model")
    scanvi_ref = scvi.model.SCANVI.from_scvi_model(
        scvi_reference_model,
        unlabeled_category=par["unknown_celltype"],
        labels_key=scvi_reference_model.adata_manager._registry["setup_args"][
            "labels_key"
        ],
    )

    reference_plan_kwargs = {
        "lr": par["learning_rate"],
        "reduce_lr_on_plateau": par["reduce_lr_on_plateau"],
        "lr_patience": par["lr_patience"],
        "lr_factor": par["lr_factor"],
    }

    logger.info("Training scANVI model on reference data with celltype labels")

    scanvi_ref.train(
        train_size=par["train_size"],
        max_epochs=par["max_epochs"],
        early_stopping=par["early_stopping"],
        early_stopping_patience=par["early_stopping_patience"],
        plan_kwargs=reference_plan_kwargs,
        check_val_every_n_epoch=1,
        accelerator="auto",
    )

    scanvi_ref.save(par["output"], overwrite=True)


if __name__ == "__main__":
    main()
