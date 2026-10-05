import sys

import mudata as mu
import numpy as np
import pandas as pd

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "rna",
    "obs_cluster": "geneformer_leiden_0.5",
    "obs_reference": None,
    "disease_clusters": ["7"],
    "disease_reference_clusters": None,
    "healthy_clusters": ["1"],
    "healthy_reference_clusters": None,
    "min_cells_per_group": 1,
    "output": "output.h5mu",
    "obs_output": "perturbation_group",
    "uns_output": "perturbation_label_cells",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402
from compress_h5mu import write_h5ad_to_h5mu_with_compression  # noqa: E402

logger = setup_logger()

LABELS = ["disease", "healthy", "unused"]


def check_column(obs, column):
    if column not in obs.columns:
        raise ValueError(
            f"Column '{column}' not found in .obs. Available: {list(obs.columns)}"
        )


def rule_mask(obs, clusters, reference_clusters):
    mask = obs[par["obs_cluster"]].astype(str).isin([str(c) for c in clusters])
    if reference_clusters:
        mask &= (
            obs[par["obs_reference"]]
            .astype(str)
            .isin([str(c) for c in reference_clusters])
        )
    return mask.to_numpy()


def main():
    logger.info("Reading modality '%s' from %s", par["modality"], par["input"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])
    obs = adata.obs

    check_column(obs, par["obs_cluster"])
    uses_reference = (
        par["disease_reference_clusters"] or par["healthy_reference_clusters"]
    )
    if uses_reference:
        if not par["obs_reference"]:
            raise ValueError(
                "Reference clusters were given, but --obs_reference was not."
            )
        check_column(obs, par["obs_reference"])

    disease = rule_mask(obs, par["disease_clusters"], par["disease_reference_clusters"])
    healthy = rule_mask(obs, par["healthy_clusters"], par["healthy_reference_clusters"])

    overlap = int((disease & healthy).sum())
    if overlap:
        raise ValueError(
            f"The disease and healthy rules both match {overlap} cells; "
            "the rules must be disjoint."
        )

    labels = np.full(adata.n_obs, "unused", dtype=object)
    labels[disease] = "disease"
    labels[healthy] = "healthy"
    adata.obs[par["obs_output"]] = pd.Categorical(labels, categories=LABELS)

    counts = {label: int((labels == label).sum()) for label in LABELS}
    logger.info("Labelled %i cells: %s", adata.n_obs, counts)
    for group in ["disease", "healthy"]:
        if counts[group] < par["min_cells_per_group"]:
            raise ValueError(
                f"The {group} rule matched {counts[group]} cells, "
                f"need at least {par['min_cells_per_group']}."
            )

    adata.uns[par["uns_output"]] = {
        "obs_cluster": par["obs_cluster"],
        "obs_reference": par["obs_reference"] or "",
        "disease_clusters": [str(c) for c in par["disease_clusters"]],
        "disease_reference_clusters": [
            str(c) for c in (par["disease_reference_clusters"] or [])
        ],
        "healthy_clusters": [str(c) for c in par["healthy_clusters"]],
        "healthy_reference_clusters": [
            str(c) for c in (par["healthy_reference_clusters"] or [])
        ],
        "n_cells": counts,
    }

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
