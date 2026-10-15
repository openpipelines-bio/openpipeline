import mudata
import sys
from harmonypy import run_harmony

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "foo.h5mu",
    "output_compression": "gzip",
    "modality": "rna",
    "obsm_input": "X_pca",
    "obsm_output": "X_pca_harmonypy",
    "obs_covariates": ["harmony_integration_leiden_1.0"],
    "theta": [2],
    "lambda": None,
    "alpha": 0.2,
    "sigma": [0.1],
    "nclust": None,
    "tau": 0,
    "block_size": 0.05,
    "max_iter_harmony": 10,
    "max_iter_kmeans": 4,
    "epsilon_cluster": 0.001,
    "epsilon_harmony": 0.01,
    "batch_prop_cutoff": 0.00001,
    "random_state": 0,
}
meta = {"resources_dir": "src/utils", "cpus": None}
### VIASH END


sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def unlist(value):
    # Pass a single value as a scalar, multiple values as a list
    if value is not None and len(value) == 1:
        return value[0]
    return value


def main():
    logger.info("Reading %s from %s", par["modality"], par["input"])
    mod = mudata.read_h5ad(par["input"].strip(), mod=par["modality"])
    logger.info("Selecting embedding %s", par["obsm_input"])
    pca_embedding = mod.obsm[par["obsm_input"]]
    metadata = mod.obs
    logger.info(
        "Available covariates to regress out: %s", ", ".join(metadata.columns.to_list())
    )
    missing_covariates = set(par["obs_covariates"]) - set(metadata.columns)
    if missing_covariates:
        raise ValueError(
            f"Covariate column(s) {', '.join(sorted(missing_covariates))} "
            f"not found in .obs of modality '{par['modality']}'."
        )
    for covariate in par["obs_covariates"]:
        n_missing = metadata[covariate].isna().sum()
        if n_missing:
            raise ValueError(
                f"Covariate column '{covariate}' contains {n_missing} missing value(s). "
                "Harmony requires every observation to have a value for each covariate."
            )

    logger.info(
        "Running harmonypy to regress out the following variables: %s",
        ", ".join(par["obs_covariates"]),
    )
    ho = run_harmony(
        pca_embedding,
        metadata,
        par["obs_covariates"],
        theta=unlist(par["theta"]),
        lamb=unlist(par["lambda"]),
        sigma=unlist(par["sigma"]),
        nclust=par["nclust"],
        tau=par["tau"],
        block_size=par["block_size"],
        max_iter_harmony=par["max_iter_harmony"],
        max_iter_kmeans=par["max_iter_kmeans"],
        epsilon_cluster=par["epsilon_cluster"],
        epsilon_harmony=par["epsilon_harmony"],
        alpha=par["alpha"],
        batch_prop_cutoff=par["batch_prop_cutoff"],
        random_state=par["random_state"],
        ncores=meta["cpus"] or 0,
    )
    logger.info("Adding output to slot %s", par["obsm_output"])
    # Since harmonypy 2.0, Z_corr is cells x PCs
    mod.obsm[par["obsm_output"]] = ho.Z_corr
    logger.info("Writing to %s", par["output"])
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], mod, par["output_compression"]
    )
    logger.info("Finished")


if __name__ == "__main__":
    main()
