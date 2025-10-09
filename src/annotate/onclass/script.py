import sys
import mudata as mu
import numpy as np
from OnClass.OnClassModel import OnClassModel
import obonet
from typing import Dict, List, Tuple


## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "model": None,
    # "reference": None,
    # "model": "resources_test/annotation_test_data/onclass_model/example_file_model",
    "reference_obs_target": "cell_ontology_class",
    "input_layer": None,
    "reference_layer": None,
    "max_iter": 100,
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "cl_nlp_emb_file": "resources_test/annotation_test_data/ontology/cl.ontology.nlp.emb",
    "cl_ontology_file": "resources_test/annotation_test_data/ontology/cl.ontology",
    "cl_obo_file": "resources_test/annotation_test_data/ontology/cl.obo",
    "output_compression": "gzip",
    "input_var_gene_names": "gene_symbol",
    "input_reference_gene_overlap": 100,
    "reference_var_input": None,
    "reference_var_gene_names": None,
    "unkown_celltype": "Unknown",
}
meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from cross_check_genes import cross_check_genes
from set_var_index import set_var_index
from subset_vars import subset_vars

logger = setup_logger()


def map_celltype_to_ontology_id(
    cl_obo_file: str,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Map cell type names to ontology IDs and vice versa.

    Parameters
    ----------
    cl_obo_file : str
        Path to the cell ontology file.

    Returns
    -------
    Tuple[Dict[str, str], Dict[str, str]]
        A tuple of two dictionaries. The first dictionary maps cell ontology IDs to cell type names.
        The second dictionary maps cell type names to cell ontology IDs.
    """
    graph = obonet.read_obo(cl_obo_file)
    cl_id_to_name = {id_: data.get("name") for id_, data in graph.nodes(data=True)}
    cl_id_to_name = {k: v for k, v in cl_id_to_name.items() if v is not None}
    name_to_cl_id = {v: k for k, v in cl_id_to_name.items()}
    return cl_id_to_name, name_to_cl_id


def cell_type_prediction(
    model: OnClassModel,
    input_matrix: np.array,
    input_features: List[str],
    id_to_name: dict,
) -> Tuple[List[str], List[float]]:
    """
    Predict cell types for input data and save results to Anndata obj.

    Parameters
    ----------
    model : OnClassModel
        The OnClass model.
    input_matrix : np.array
        The input data matrix.
    input_modality : ad.AnnData
        The input data Anndata object.
    id_to_name : dict
        Dictionary mapping cell ontology IDs to cell type names.

    Returns
    -------
    predictions: List[str]
        The predicted cell types.
    probabilities: List[float]
        Probabilities of the predicted cell types.
    """
    corr_test_feature = model.ProcessTestFeature(
        test_feature=input_matrix,
        test_genes=input_features,
        log_transform=False,
    )
    onclass_pred = model.Predict(
        corr_test_feature, use_normalize=False, refine=True, unseen_ratio=-1.0
    )
    pred_label = [model.i2co[ind] for ind in onclass_pred[2]]
    pred_cell_type_label = [id_to_name[id] for id in pred_label]
    prob_cell_type_label = np.max(onclass_pred[1], axis=1) / onclass_pred[1].sum(1)

    return pred_cell_type_label, prob_cell_type_label


def main():
    if (not par["model"] and not par["reference"]) or (
        par["model"] and par["reference"]
    ):
        raise ValueError(
            "Make sure to provide either 'model' or 'reference', but not both."
        )

    logger.info("Reading input data")
    input_mudata = mu.read_h5mu(par["input"])
    input_adata = input_mudata.mod[par["modality"]]
    input_modality = input_adata.copy()

    # Set var names to the desired gene name format (gene symbol, ensembl id, etc.)
    input_modality = set_var_index(
        input_modality, par["input_var_gene_names"], par["sanitize_gene_names"]
    )
    input_matrix = (
        input_modality.layers[par["input_layer"]]
        if par["input_layer"]
        else input_modality.X
    )
    # Onclass needs dense matrix format
    input_matrix = input_matrix.toarray()

    id_to_name, name_to_id = map_celltype_to_ontology_id(par["cl_obo_file"])

    if par["model"]:
        logger.info("Predicting cell types using pre-trained model")
        model = OnClassModel(
            cell_type_nlp_emb_file=par["cl_nlp_emb_file"],
            cell_type_network_file=par["cl_ontology_file"],
        )

        model.BuildModel(use_pretrain=par["model"], ngene=None)
        cross_check_genes(
            model.genes, input_modality.var.index, par["input_reference_gene_overlap"]
        )

    elif par["reference"]:
        logger.info("Reading reference data")
        model = OnClassModel(
            cell_type_nlp_emb_file=par["cl_nlp_emb_file"],
            cell_type_network_file=par["cl_ontology_file"],
        )

        reference_mudata = mu.read_h5mu(par["reference"])
        reference_modality = reference_mudata.mod[par["modality"]].copy()
        reference_modality = set_var_index(
            reference_modality,
            par["reference_var_gene_names"],
            par["sanitize_gene_names"],
        )

        # subset to HVG if required
        if par["reference_var_input"]:
            reference_modality = subset_vars(
                reference_modality, par["reference_var_input"]
            )

        cross_check_genes(
            input_modality.var.index,
            reference_modality.var.index,
            par["input_reference_gene_overlap"],
        )

        reference_matrix = (
            reference_modality.layers[par["reference_layer"]]
            if par["reference_layer"]
            else reference_modality.X
        )
        # Onclass needs dense matrix format
        reference_matrix = reference_matrix.toarray()

        logger.info("Training a model from reference...")

        labels = reference_modality.obs[par["reference_obs_target"]].tolist()
        labels_cl = [
            name_to_id[label] if label in name_to_id else par["unknown_celltype"]
            for label in labels
        ]

        _ = model.EmbedCellTypes(labels_cl)
        corr_train_feature, _, corr_train_genes, _ = model.ProcessTrainFeature(
            train_feature=reference_matrix,
            train_label=labels_cl,
            train_genes=reference_modality.var.index,
            test_feature=input_matrix,
            test_genes=input_modality.var.index,
            log_transform=False,
        )
        model.BuildModel(ngene=len(corr_train_genes))
        model.Train(corr_train_feature, labels_cl, max_iter=par["max_iter"])

    logger.info("Predicting cell types")
    predictions, probabilities = cell_type_prediction(
        model, input_matrix, input_modality.var.index, id_to_name
    )

    logger.info("Writing output data")
    input_adata.obs[par["output_obs_predictions"]] = predictions
    input_adata.obs[par["output_obs_probability"]] = probabilities
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == "__main__":
    main()
