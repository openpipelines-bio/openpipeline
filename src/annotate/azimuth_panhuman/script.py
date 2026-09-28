import sys

import mudata as mu
import pandas as pd
from scipy.sparse import csr_matrix

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "input_layer": None,
    "input_var_gene_names": "gene_symbol",
    "sanitize_ensembl_ids": False,
    "input_reference_gene_overlap": 100,
    "model_version": "v1",
    "annotation_pipeline": "supervised",
    "eval_batch_size": 8192,
    "normalization_override": False,
    "norm_check_batch_size": 100,
    "output_mode": "minimal",
    "refine_labels": True,
    "map_to_cl": None,
    "include_cl_id": False,
    "extract_embeddings": True,
    "umap_embeddings": True,
    "umap_n_neighbors": 30,
    "umap_n_components": 2,
    "umap_metric": "cosine",
    "umap_min_dist": 0.3,
    "umap_lr": 1.0,
    "umap_seed": 42,
    "umap_spread": 1.0,
    "umap_init": "spectral",
    "umap_verbose": False,
    "output_obsm_embedding": "X_azimuth",
    "output_obsm_umap": "X_azimuth_umap",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from cross_check_genes import cross_check_genes
from set_var_index import set_var_index
from setup_logger import setup_logger

logger = setup_logger()

import panhumanpy as ph
from panhumanpy.ANNotate_tools import InferenceTools, check_normalization


def main(par):
    logger.info("Reading input data")
    input_mudata = mu.read_h5mu(par["input"])
    input_adata = input_mudata.mod[par["modality"]]

    query_adata = set_var_index(
        input_adata.copy(), par["input_var_gene_names"], par["sanitize_ensembl_ids"]
    )

    count_matrix = (
        query_adata.layers[par["input_layer"]] if par["input_layer"] else query_adata.X
    )
    X_query = csr_matrix(count_matrix)
    query_features = query_adata.var.index.astype(str).tolist()

    # panhumanpy itself never raises on this: it only heuristically guesses
    # whether the data is already normalized and silently proceeds either
    # way. Fail loudly instead, matching celltypist/singler's convention,
    # unless the user explicitly asked to skip the check.
    if not par["normalization_override"] and check_normalization(
        X_query, par["normalization_override"], par["norm_check_batch_size"]
    ):
        raise ValueError(
            "Invalid expression matrix: detected non-integer values in "
            "--input_layer (or .X if not set), suggesting the data is "
            "already normalized. Azimuth expects raw counts and performs "
            "its own normalization internally. Pass --normalization_override "
            "if you are certain this is a false positive."
        )

    # Only reads the (package-bundled) reference gene panel, not the
    # downloaded neural network weights, so this stays cheap even though
    # annotate_core() below loads the full model a second time.
    logger.info("Checking gene overlap with the Azimuth reference gene panel")
    feature_panel = InferenceTools(
        annotation_pipeline=par["annotation_pipeline"],
        model_version=par["model_version"],
    ).load_inference_feature_panel()
    cross_check_genes(
        query_features,
        feature_panel,
        min_gene_overlap=par["input_reference_gene_overlap"],
    )

    cells_meta = pd.DataFrame(index=query_adata.obs_names)

    logger.info("Running Azimuth annotation")
    core_outputs = ph.annotate_core(
        X_query,
        query_features,
        cells_meta,
        annotation_pipeline=par["annotation_pipeline"],
        eval_batch_size=par["eval_batch_size"],
        normalization_override=par["normalization_override"],
        norm_check_batch_size=par["norm_check_batch_size"],
        output_mode=par["output_mode"],
        refine_labels=par["refine_labels"],
        map_to_cl=par["map_to_cl"],
        include_cl_id=par["include_cl_id"],
        extract_embeddings=par["extract_embeddings"],
        umap_embeddings=par["umap_embeddings"],
        n_neighbors=par["umap_n_neighbors"],
        n_components=par["umap_n_components"],
        metric=par["umap_metric"],
        min_dist=par["umap_min_dist"],
        umap_lr=par["umap_lr"],
        umap_seed=par["umap_seed"],
        spread=par["umap_spread"],
        verbose=par["umap_verbose"],
        init=par["umap_init"],
        model_version=par["model_version"],
    )

    cells_meta_out = core_outputs["cells_meta"]
    embeddings_dict = core_outputs["embeddings_dict"]
    umap_dict = core_outputs["umap_dict"]

    logger.info("Writing annotations to output object")
    for col in cells_meta_out.columns:
        input_adata.obs[col] = cells_meta_out[col].values

    if par["extract_embeddings"] and "azimuth_embed" in embeddings_dict:
        input_adata.obsm[par["output_obsm_embedding"]] = embeddings_dict[
            "azimuth_embed"
        ]

    if par["umap_embeddings"] and "azimuth_umap" in umap_dict:
        input_adata.obsm[par["output_obsm_umap"]] = umap_dict["azimuth_umap"]

    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main(par)
