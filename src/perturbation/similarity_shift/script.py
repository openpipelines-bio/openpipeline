import sys

import mudata as mu
import numpy as np
import pandas as pd

## VIASH START
par = {
    "input": "virtual_cells.h5mu",
    "original": "original_cells.h5mu",
    "modality": "rna",
    "obsm_input": "X_geneformer",
    "obs_cell_id": "perturbation_cell_id",
    "obs_gene_id": "perturbation_gene_id",
    "obs_gene_name": "perturbation_gene_name",
    "uns_centroids": "perturbation_centroids",
    "healthy_label": "healthy",
    "disease_label": "disease",
    "output": "similarity_shift.csv",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger  # noqa: E402

logger = setup_logger()

EPS = 1e-6


def cosine_to(matrix, centroid):
    # row-wise cosine similarity of every row of matrix to one vector
    norms = np.linalg.norm(matrix, axis=1) * np.linalg.norm(centroid)
    return (matrix @ centroid) / np.clip(norms, EPS, None)


def get_embedding(adata, path):
    if par["obsm_input"] not in adata.obsm:
        raise ValueError(
            f"'{par['obsm_input']}' not found in .obsm of {path}. "
            f"Available: {list(adata.obsm)}"
        )
    return np.asarray(adata.obsm[par["obsm_input"]], dtype=np.float64)


def get_centroids(original):
    if par["uns_centroids"] not in original.uns:
        raise ValueError(
            f"'{par['uns_centroids']}' not found in .uns of {par['original']}. "
            "Compute the centroids with perturbation/compute_centroids first."
        )
    centroids = original.uns[par["uns_centroids"]]
    missing = [
        label
        for label in [par["healthy_label"], par["disease_label"]]
        if label not in centroids.index
    ]
    if missing:
        raise ValueError(
            f"No centroid for {missing} in .uns['{par['uns_centroids']}']. "
            f"Available: {list(centroids.index)}"
        )
    return (
        centroids.loc[par["healthy_label"]].to_numpy(dtype=np.float64),
        centroids.loc[par["disease_label"]].to_numpy(dtype=np.float64),
    )


def main():
    virtual = mu.read_h5ad(par["input"], mod=par["modality"])
    original = mu.read_h5ad(par["original"], mod=par["modality"])

    for column in [par["obs_cell_id"], par["obs_gene_id"], par["obs_gene_name"]]:
        if column not in virtual.obs.columns:
            raise ValueError(f"Column '{column}' not found in .obs of {par['input']}")

    perturbed = get_embedding(virtual, par["input"])
    healthy, disease = get_centroids(original)
    if perturbed.shape[1] != healthy.shape[0]:
        raise ValueError(
            f"The embeddings have {perturbed.shape[1]} dimensions, "
            f"the centroids {healthy.shape[0]}."
        )

    cell_ids = virtual.obs[par["obs_cell_id"]].astype(str)
    unknown = sorted(set(cell_ids) - set(original.obs_names))
    if unknown:
        raise ValueError(
            f"{len(unknown)} cells of {par['input']} are not in {par['original']}, "
            f"e.g. {unknown[:5]}"
        )
    original_embedding = get_embedding(original, par["original"])[
        original.obs_names.get_indexer(cell_ids)
    ]

    result = pd.DataFrame(
        {
            "cell_id": cell_ids.to_numpy(),
            "gene_id": virtual.obs[par["obs_gene_id"]].astype(str).to_numpy(),
            "gene_name": virtual.obs[par["obs_gene_name"]].astype(str).to_numpy(),
            "cos_original_healthy": cosine_to(original_embedding, healthy),
            "cos_original_disease": cosine_to(original_embedding, disease),
            "cos_perturbed_healthy": cosine_to(perturbed, healthy),
            "cos_perturbed_disease": cosine_to(perturbed, disease),
        }
    )
    result["shift_healthy"] = (
        result["cos_perturbed_healthy"] - result["cos_original_healthy"]
    )
    result["shift_disease"] = (
        result["cos_perturbed_disease"] - result["cos_original_disease"]
    )

    result.to_csv(par["output"], index=False)
    logger.info(
        "Scored %i virtual cells of %i cells, written to %s",
        result.shape[0],
        cell_ids.nunique(),
        par["output"],
    )


if __name__ == "__main__":
    main()
