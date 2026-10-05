import sys

import mudata as mu
import numpy as np
import pandas as pd

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "rna",
    "obsm_input": "X_geneformer",
    "obs_group": "perturbation_group",
    "obs_filter": None,
    "groups": ["disease", "healthy"],
    "method": "median",
    "min_cells": 1,
    "output": "output.h5mu",
    "uns_output": "perturbation_centroids",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402
from compress_h5mu import write_h5ad_to_h5mu_with_compression  # noqa: E402

logger = setup_logger()


def main():
    logger.info("Reading modality '%s' from %s", par["modality"], par["input"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])

    if par["obsm_input"] not in adata.obsm:
        raise ValueError(
            f"'{par['obsm_input']}' not found in .obsm. Available: {list(adata.obsm)}"
        )
    for column in [par["obs_group"], par["obs_filter"]]:
        if column and column not in adata.obs.columns:
            raise ValueError(
                f"Column '{column}' not found in .obs. Available: {list(adata.obs.columns)}"
            )

    embedding = np.asarray(adata.obsm[par["obsm_input"]], dtype=np.float64)
    groups = adata.obs[par["obs_group"]].astype(str).to_numpy()
    keep = np.ones(adata.n_obs, dtype=bool)
    if par["obs_filter"]:
        keep = adata.obs[par["obs_filter"]].astype(bool).to_numpy()
        logger.info("%i of %i cells pass --obs_filter", keep.sum(), adata.n_obs)

    wanted = par["groups"] or sorted(set(groups[keep]))
    summarize = np.median if par["method"] == "median" else np.mean

    centroids = {}
    for group in wanted:
        cells = keep & (groups == group)
        if cells.sum() < par["min_cells"]:
            raise ValueError(
                f"Group '{group}' has {cells.sum()} cells, need at least {par['min_cells']}."
            )
        centroids[group] = summarize(embedding[cells], axis=0)
        logger.info(
            "%s centroid of '%s' over %i cells", par["method"], group, cells.sum()
        )

    adata.uns[par["uns_output"]] = pd.DataFrame.from_dict(
        centroids,
        orient="index",
        columns=[str(i) for i in range(embedding.shape[1])],
    )

    logger.info("Writing output to %s", par["output"])
    write_h5ad_to_h5mu_with_compression(
        output_file=par["output"],
        h5mu=par["input"],
        modality_name=par["modality"],
        modality_data=adata,
        output_compression=par["output_compression"],
    )


if __name__ == "__main__":
    main()
