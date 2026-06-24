from __future__ import annotations
import sys
import numpy as np

## VIASH START
par = {
    "input": "test_data/phate_input.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca",
    "output": "test_data/phate_output.h5mu",
    "obsm_output": "X_phate",
    "n_components": 2,
    "knn": 5,
    "decay": 40,
    "t": "auto",
    "gamma": 1.0,
    "random_state": 0,
    "output_compression": None,
}
meta = {
    "resources_dir": ".",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def main():
    import mudata as mu
    import phate

    logger.info("Reading '%s', modality '%s'.", par["input"], par["modality"])
    data = mu.read_h5ad(par["input"], mod=par["modality"])

    obsm_key = par["obsm_input"]
    if obsm_key not in data.obsm:
        raise ValueError(
            f"'{obsm_key}' not found in .mod['{par['modality']}'].obsm. "
            f"Available keys: {list(data.obsm.keys())}"
        )

    X = np.array(data.obsm[obsm_key])
    logger.info(
        "Running PHATE on '%s' (shape %s): n_components=%d, knn=%d, decay=%d, t=%s.",
        obsm_key,
        X.shape,
        par["n_components"],
        par["knn"],
        par["decay"],
        par["t"],
    )

    # Parse t: keep as "auto" string or convert to int
    t = par["t"]
    if t != "auto":
        try:
            t = int(t)
        except ValueError:
            raise ValueError(f"--t must be 'auto' or a positive integer, got '{t}'.")

    n_jobs = meta["cpus"] if meta.get("cpus") else 1

    phate_op = phate.PHATE(
        n_components=par["n_components"],
        knn=par["knn"],
        decay=par["decay"],
        t=t,
        gamma=par["gamma"],
        random_state=par["random_state"],
        n_jobs=n_jobs,
        verbose=False,
    )

    X_phate = phate_op.fit_transform(X)

    data.obsm[par["obsm_output"]] = X_phate
    logger.info(
        "Stored PHATE embedding in .obsm['%s'] (shape %s).",
        par["obsm_output"],
        X_phate.shape,
    )

    logger.info(
        "Writing output to '%s' with compression '%s'.",
        par["output"],
        par["output_compression"],
    )
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], data, par["output_compression"]
    )
    logger.info("Finished.")


if __name__ == "__main__":
    main()
