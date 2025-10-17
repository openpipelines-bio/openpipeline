import sys
import celltypist
import mudata as mu
import anndata as ad
import pandas as pd
import scanpy as sc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "model": None,
    "input_layer": "log_normalized",
    "reference_layer": "log_normalized",
    "input_reference_gene_overlap": 100,
    "reference_obs_target": "cell_ontology_class",
    "reference_var_input": None,
    "feature_selection": True,
    "majority_voting": True,
    "output_compression": "gzip",
    "input_var_gene_names": None,
    "reference_var_gene_names": "ensemblid",
    "C": 1.0,
    "max_iter": 1000,
    "use_SGD": False,
    "min_prop": 0,
    "output_obs_predictions": "celltypist_pred",
    "output_obs_probability": "celltypist_probability",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from cross_check_genes import cross_check_genes
from set_var_index import set_var_index
from subset_vars import subset_vars

logger = setup_logger()


def setup_anndata(
    adata: ad.AnnData,
    layer: str | None = None,
    gene_names: str | None = None,
    sanitize_gene_names: bool = True,
    var_input: str | None = None,
) -> ad.AnnData:
    """Creates an AnnData object in the expected format for CellTypist,
    with lognormalized data (with a target sum of 10000) in the .X slot.

    Parameters
    ----------
    adata
        AnnData object.
    layer
        Layer in AnnData object to lognormalize.
    gene_names
        .obs field with the gene names to be used
    var_input
        .var field with a boolean array of the genes to be used (e.g. highly variable genes)
    Returns
    -------
    AnnData object in CellTypist format.
    """

    adata = set_var_index(adata, gene_names, sanitize_gene_names)

    if var_input:
        adata = subset_vars(adata, var_input)

    raw_counts = adata.layers[layer].copy() if layer else adata.X.copy()

    input_modality = ad.AnnData(X=raw_counts, var=pd.DataFrame(index=adata.var.index))

    sc.pp.normalize_total(input_modality, target_sum=10000)
    sc.pp.log1p(input_modality)

    return input_modality


def main(par):
    if (not par["model"] and not par["reference"]) or (
        par["model"] and par["reference"]
    ):
        raise ValueError(
            "Make sure to provide either 'model' or 'reference', but not both."
        )

    input_mudata = mu.read_h5mu(par["input"])
    input_adata = input_mudata.mod[par["modality"]]
    input_modality = input_adata.copy()

    # Provide correct format of query data for celltypist annotation
    input_modality = setup_anndata(
        input_modality,
        par["input_layer"],
        par["input_var_gene_names"],
        par["sanitize_gene_names"],
    )

    if par["model"]:
        logger.info("Loading CellTypist model")
        model = celltypist.models.Model.load(par["model"])
        cross_check_genes(
            input_modality.var.index,
            model.features,
            min_gene_overlap=par["input_reference_gene_overlap"],
        )

    elif par["reference"]:
        reference_adata = mu.read_h5mu(par["reference"]).mod[par["modality"]]
        reference_modality = reference_adata.copy()

        # Provide correct format of query data for celltypist annotation
        reference_modality = setup_anndata(
            reference_modality,
            par["reference_layer"],
            par["reference_var_gene_names"],
            par["sanitize_gene_names"],
            par["reference_var_input"],
        )

        # Ensure enough overlap between genes in query and reference
        cross_check_genes(
            input_modality.var.index,
            reference_modality.var.index,
            min_gene_overlap=par["input_reference_gene_overlap"],
        )

        logger.info("Training CellTypist model on reference")
        model = celltypist.train(
            reference_modality.X,
            labels=reference_adata.obs[par["reference_obs_target"]],
            genes=reference_modality.var.index,
            C=par["C"],
            max_iter=par["max_iter"],
            use_SGD=par["use_SGD"],
            feature_selection=par["feature_selection"],
            check_expression=True,
        )

    logger.info("Predicting CellTypist annotations")
    predictions = celltypist.annotate(
        input_modality, model, majority_voting=par["majority_voting"]
    )

    input_adata.obs[par["output_obs_predictions"]] = predictions.predicted_labels[
        "predicted_labels"
    ].values
    input_adata.obs[par["output_obs_probability"]] = predictions.probability_matrix.max(
        axis=1
    ).values

    # copy observations back to input data (with full set of features)
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main(par)
