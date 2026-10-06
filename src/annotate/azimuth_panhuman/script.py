import sys
import tempfile
from pathlib import Path

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
    "input_reference_gene_overlap": 100,
    "model": None,
    "model_version": "v1",
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
    "output_obs_predictions": "azimuth_pred",
    "output_obs_probability": "azimuth_probability",
    "output_obs_predictions_broad": "azimuth_broad",
    "output_obs_predictions_medium": "azimuth_medium",
    "output_obs_predictions_fine": "azimuth_fine",
    "output_obsm_embedding": "X_azimuth",
    "output_obsm_umap": "X_azimuth_umap",
    "output_compression": None,
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from cross_check_genes import cross_check_genes
from is_lognormalized import is_lognormalized
from set_var_index import set_var_index
from setup_logger import setup_logger

logger = setup_logger()

import panhumanpy as ph
import panhumanpy.ANNotate_tools as ANNotate_tools
import tensorflow as tf
from panhumanpy.ANNotate_tools import InferenceTools, check_normalization

# Currently the only annotation pipeline implemented by panhumanpy.
ANNOTATION_PIPELINE = "supervised"


def stage_model(model_path, model_version):
    """Place the model where panhumanpy's loader looks for it, so it is not downloaded.

    panhumanpy cannot be given a model path directly: it loads
    <CACHE_DIR>/<model_version>/inference_model/inference_model.keras and only
    downloads the weights when that file is missing. Point CACHE_DIR at a
    temporary directory containing a symlink to the provided model.
    """
    cache_dir = Path(tempfile.mkdtemp())
    model_dir = cache_dir / model_version / "inference_model"
    model_dir.mkdir(parents=True)
    (model_dir / "inference_model.keras").symlink_to(Path(model_path).resolve())
    ANNotate_tools.CACHE_DIR = cache_dir


def main(par):
    gpu_devices = tf.config.list_physical_devices("GPU")
    logger.info(
        "GPU devices visible to TensorFlow: %s",
        gpu_devices if gpu_devices else "none, running on CPU",
    )

    if par["model"]:
        logger.info("Using provided model %s", par["model"])
        stage_model(par["model"], par["model_version"])

    logger.info("Reading input data")
    input_mudata = mu.read_h5mu(par["input"])
    input_adata = input_mudata.mod[par["modality"]]

    # Azimuth expects gene symbols, not Ensembl IDs (see --input_var_gene_names),
    # so there is nothing to sanitize.
    query_adata = set_var_index(
        input_adata.copy(), par["input_var_gene_names"], sanitize_ensembl_ids=False
    )

    count_matrix = (
        query_adata.layers[par["input_layer"]] if par["input_layer"] else query_adata.X
    )
    X_query = csr_matrix(count_matrix)
    query_features = query_adata.var.index.astype(str).tolist()

    if par["normalization_override"]:
        # Azimuth's internal normalization is skipped when this flag is set,
        # so verify the data was already normalized the way Azimuth expects
        # (log1p to a target sum of 10000 counts per cell)
        if not is_lognormalized(X_query, target_sum=10000):
            raise ValueError(
                "Invalid expression matrix: --normalization_override was "
                "set, but --input_layer (or .X if not set) does not look "
                "like it was log1p-normalized to a target sum of 10000 "
                "counts per cell, which is what Azimuth expects when its "
                "internal normalization is skipped."
            )
    else:
        # panhumanpy itself never raises on this: it only heuristically
        # guesses whether the data is already normalized and silently
        # proceeds either way. Fail loudly instead.
        if check_normalization(
            X_query, par["normalization_override"], par["norm_check_batch_size"]
        ):
            raise ValueError(
                "Invalid expression matrix: detected non-integer values in "
                "--input_layer (or .X if not set), suggesting the data is "
                "already normalized. Azimuth expects raw counts and performs "
                "its own normalization internally. Pass "
                "--normalization_override if the data is already "
                "log1p-normalized to a target sum of 10000 counts per cell."
            )

    # Only reads the (package-bundled) reference gene panel, not the
    # downloaded neural network weights, so this stays cheap even though
    # annotate_core() below loads the full model a second time.
    logger.info("Checking gene overlap with the Azimuth reference gene panel")
    feature_panel = InferenceTools(
        annotation_pipeline=ANNOTATION_PIPELINE,
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
        annotation_pipeline=ANNOTATION_PIPELINE,
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
    # Map the Azimuth output columns to their parametrized .obs column names.
    # Only these columns are copied to the input data, so no other existing
    # .obs columns can be overwritten.
    output_obs_columns = {
        "final_level_labels": par["output_obs_predictions"],
        "final_level_confidence": par["output_obs_probability"],
    }
    if par["refine_labels"]:
        output_obs_columns.update(
            {
                "azimuth_broad": par["output_obs_predictions_broad"],
                "azimuth_medium": par["output_obs_predictions_medium"],
                "azimuth_fine": par["output_obs_predictions_fine"],
            }
        )
    for azimuth_col, output_col in output_obs_columns.items():
        input_adata.obs[output_col] = cells_meta_out[azimuth_col].values

    if par["extract_embeddings"] and "azimuth_embed" in embeddings_dict:
        input_adata.obsm[par["output_obsm_embedding"]] = embeddings_dict[
            "azimuth_embed"
        ]

    if par["umap_embeddings"] and "azimuth_umap" in umap_dict:
        input_adata.obsm[par["output_obsm_umap"]] = umap_dict["azimuth_umap"]

    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main(par)
