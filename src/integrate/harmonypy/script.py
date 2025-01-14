import mudata
import sys
from harmonypy import run_harmony

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "foo.h5mu",
    "compression": "gzip",
    "modality": "rna",
    "obsm_input": "X_pca",
    "obsm_output": "X_pca_harmonypy",
    "theta": 2,
    "obs_covariates": ["batch"],
}
### VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def main():
    logger.info("Reading %s from %s", par["modality"], par["input"])
    mod = mudata.read_h5ad(par["input"].strip(), mod=par["modality"])
    logger.info("Selecting embedding %s", par["obsm_input"])
    pca_embedding = mod.obsm[par["obsm_input"]]
    metadata = mod.obs
    logger.info(
        "Available covariates to regress out: ", ", ".join(metadata.columns.to_list())
    )
    logger.info(
        "Running harmonypy with theta %s to regress out the following variables: ",
        par["theta"],
        ", ".join(par["obs_covariates"]),
    )
    ho = run_harmony(pca_embedding, metadata, par["obs_covariates"], theta=par["theta"])
    logger.info("Adding output to slot %s", par["obsm_output"])
    mod.obsm[par["obsm_output"]] = ho.Z_corr.T
    logger.info("Writing to %s")
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], mod, par["output_compression"]
    )
    logger.info("Finished")


if __name__ == "__main__":
    main()
