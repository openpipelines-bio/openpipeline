import sys

import cellmapper
import mudata as mu
import numpy as np

## VIASH START
par = {
    "input": "query.h5mu",
    "modality": "rna",
    "input_obsm_features": "X_pca",
    "reference": "reference.h5mu",
    "reference_obsm_features": "X_pca",
    "reference_obs_targets": ["cell_type"],
    "output": "output.h5mu",
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "output_compression": "gzip",
    "n_neighbors": 15,
    "kernel_method": None,
}
meta = {"resources_dir": "src/labels_transfer/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from helper import check_arguments, get_query_features, get_reference_features
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression


logger = setup_logger()


def validate_inputs(par, query_adata, ref_adata):
    modality = par["modality"]

    embeddings = {
        "query": (query_adata, par["input_obsm_features"]),
        "reference": (ref_adata, par["reference_obsm_features"]),
    }

    for name, (adata, key) in embeddings.items():
        if key is not None and key not in adata.obsm:
            raise ValueError(
                f"Embedding '{key}' not found in {name} modality '{modality}' obsm. "
                f"Available keys: {list(adata.obsm.keys())}"
            )

    for target in par["reference_obs_targets"]:
        if target not in ref_adata.obs:
            raise ValueError(
                f"Reference label '{target}' not found in reference modality '{modality}' obs. "
                f"Available keys: {list(ref_adata.obs.keys())}"
            )


def main(par):
    par = check_arguments(par)

    logger.info(
        "Loading query AnnData from '%s' modality '%s'", par["input"], par["modality"]
    )
    query_adata = mu.read_anndata(par["input"], mod=par["modality"])
    logger.info(
        "Loading reference AnnData from '%s' modality '%s'",
        par["reference"],
        par["modality"],
    )
    ref_adata = mu.read_anndata(par["reference"], mod=par["modality"])

    validate_inputs(par, query_adata, ref_adata)

    query_features = get_query_features(query_adata, par, logger)
    ref_features = get_reference_features(ref_adata, par, logger)

    q_dims = query_features.shape[1]
    r_dims = ref_features.shape[1]
    if q_dims != r_dims:
        raise ValueError(
            f"Embedding dimensions do not match! Query: {q_dims}, Reference: {r_dims}. "
            "Reference and query must share the same embedding space."
        )

    # CellMapper needs both embeddings to have the same key
    use_rep = "_X_cellmapper_features"
    query_adata.obsm[use_rep] = query_features
    ref_adata.obsm[use_rep] = ref_features

    logger.info("Running CellMapper")
    mapper = cellmapper.CellMapper(query_adata, ref_adata)
    mapper.map(
        obs_keys=par["reference_obs_targets"],
        use_rep=use_rep,
        n_neighbors=par["n_neighbors"],
        prediction_postfix="_pred",
        kernel_method=par["kernel_method"],
    )

    # Map CellMapper's prediction column names to the user-specified output column names
    for target, output_pred, output_prob in zip(
        par["reference_obs_targets"],
        par["output_obs_predictions"],
        par["output_obs_probability"],
    ):
        default_pred_col = f"{target}_pred"
        if default_pred_col in query_adata.obs and output_pred != default_pred_col:
            query_adata.obs[output_pred] = query_adata.obs.pop(default_pred_col)

        if output_pred not in query_adata.obs and default_pred_col in query_adata.obs:
            query_adata.obs[output_pred] = query_adata.obs[default_pred_col]

        if output_prob not in query_adata.obs:
            # CellMapper does not expose per-label probabilities in the returned obs.
            query_adata.obs[output_prob] = np.nan

    # Remove the temporary embedding
    query_adata.obsm.pop(use_rep, None)
    logger.info("Writing output to '%s'", par["output"])
    write_h5ad_to_h5mu_with_compression(
        par["output"],
        par["input"],
        par["modality"],
        query_adata,
        par["output_compression"],
    )


if __name__ == "__main__":
    sys.exit(main(par))
