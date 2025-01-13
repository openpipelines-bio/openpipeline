import sys
import mudata as mu
import re


## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "input_layer": None,
    "input_var_gene_names": None,
    "input_obs_batch": "sample_id",
    "input_obs_labels": None,
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "reference_layer": None,
    "reference_var_gene_names": "ensemblid",
    "reference_obs_batch": "sample_id",
    "reference_obs_label": "cell_ontology_class",
    "output_query": "pbmc_1k_protein_v3_filtered_feature_bc_matrix_aligned.h5mu",
    "output_reference": "TS_Blood_filtered_aligned.h5mu",
    "output_layer": "_counts",
    "output_var_gene_names": "_gene_names",
    "output_obs_batch": "_sample_id",
    "output_obs_label": "_cell_type",
    "input_reference_gene_overlap": 100,
    "unkown_celltype_label": "Unknown",
    "overwrite_existing_key": False,
}

meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from cross_check_genes import cross_check_genes

logger = setup_logger()


def copy_layer(adata, input_layer, output_layer, overwrite=False):
    if output_layer not in adata.layers.keys() and input_layer:
        logger.info(f"Copying {input_layer} layer to {output_layer}")
        adata.layers[output_layer] = adata.layers[output_layer].copy()

    elif output_layer not in adata.layers.keys() and input_layer:
        logger.info(f"Copying .X layer to {output_layer}")
        adata.layers[output_layer] = adata.X.copy()

    else:
        if not overwrite:
            raise ValueError(
                f"Layer {output_layer} already exists. Data can not be copied."
            )
        logger.warning(f"Layer {output_layer} already exists. Overwriting data.")
        adata.layers[output_layer] = (
            adata.layers[input_layer].copy() if input_layer else adata.X.copy()
        )

    return adata


def copy_obs(adata, input_obs_key, output_obs_key, overwrite=False, fill_value=None):
    if output_obs_key not in adata.obs.keys():
        logger.info(f"Copying {input_obs_key} .obs key to {output_obs_key}")
        adata.layers[output_obs_key] = (
            adata.layers[input_obs_key].copy() if input_obs_key else fill_value
        )

    else:
        if not overwrite:
            raise ValueError(
                f".obs key {output_obs_key} already exists. Data can not be copied."
            )

        logger.warning(f".obs key {output_obs_key} already exists. Overwriting data.")
        adata.layers[output_obs_key] = (
            adata.layers[input_obs_key].copy() if input_obs_key else fill_value
        )

    return adata


def copy_and_sanitize_var_gene_names(
    adata, input_var_key, output_var_key, overwrite=False
):
    if output_var_key not in adata.var.keys() and input_var_key:
        logger.info(f"Copying {input_var_key} .var key to {output_var_key}")
        adata.var[output_var_key] = [
            re.sub("\\.[0-9]+$", "", s) for s in adata.var[input_var_key]
        ]
    if output_var_key not in adata.var.keys() and input_var_key:
        logger.info(f"Copying {input_var_key} .var key to {output_var_key}")
        adata.var[output_var_key] = [
            re.sub("\\.[0-9]+$", "", s) for s in adata.var.index
        ]

    else:
        if not overwrite:
            raise ValueError(
                f".var key {output_var_key} already exists. Data can not be copied."
            )

        logger.warning(f".var key {output_var_key} already exists. Overwriting data.")
        adata.var[output_var_key] = [
            re.sub("\\.[0-9]+$", "", s) for s in adata.var[input_var_key]
        ]

    return adata


def main():
    logger.info("Reading query and reference data")
    input_mudata = mu.read_h5mu(par["input"])
    input_adata = input_mudata.mod[par["modality"]]
    input_modality = input_adata.copy()

    reference_mudata = mu.read_h5mu(par["reference"])
    reference_adata = reference_mudata.mod[par["modality"]]
    reference_modality = reference_adata.copy()

    # Aligning layers
    logger.info("Aligning layers")
    input_modality = copy_layer(
        input_modality,
        par["input_layer"],
        par["output_layer"],
        overwrite=par["overwrite_existing_key"],
    )
    reference_modality = copy_layer(
        reference_modality,
        par["reference_layer"],
        par["output_layer"],
        overwrite=par["overwrite_existing_key"],
    )

    # Aligning batch labels
    logger.info("Aligning batch labels")
    input_modality = copy_obs(
        input_modality,
        par["input_obs_batch"],
        par["output_obs_batch"],
        par["overwrite_existing_key"],
        fill_value=None,
    )
    reference_modality = copy_obs(
        reference_modality,
        par["reference_obs_batch"],
        par["output_obs_batch"],
        overwrite=par["overwrite_existing_key"],
        fill_value=None,
    )

    # Aligning cell type labels
    logger.info("Aligning cell type labels")
    input_modality = copy_obs(
        input_modality,
        par["input_obs_celltype_label"],
        par["output_obs_celltype_label"],
        par["overwrite_existing_key"],
        fill_vallue=par["unkown_celltype_label"],
    )
    reference_modality = copy_obs(
        reference_modality,
        par["reference_obs_celltype_label"],
        par["output_obs_celltype_label"],
        overwrite=par["overwrite_existing_key"],
        fill_vallue=par["unkown_celltype_label"],
    )

    # Aligning and sanitizing gene names
    logger.info("Aligning and sanitizing gene names")
    input_modality = copy_and_sanitize_var_gene_names(
        input_modality,
        par["input_var_gene_names"],
        par["output_var_gene_names"],
        overwrite=par["overwrite_existing_key"],
    )
    reference_modality = copy_and_sanitize_var_gene_names(
        reference_modality,
        par["reference_var_gene_names"],
        par["output_var_gene_names"],
        overwrite=par["overwrite_existing_key"],
    )

    # Cross check genes
    cross_check_genes(
        input_modality.var["output_var_gene_names"],
        reference_modality.var["output_var_gene_names"],
        min_gene_overlap=par["input_reference_gene_overlap"],
    )

    # Writing output data
    logger.info("Writing output data")
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output_query"], compression=par["output_compression"])

    reference_mudata.mod[par["modality"]] = reference_modality
    reference_mudata.write_h5mu(
        par["output_reference"], compression=par["output_compression"]
    )
