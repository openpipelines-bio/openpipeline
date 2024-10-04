import sys
import logging
import mudata as mu
import anndata as ad
import re
import numpy as np
from OnClass.OnClassModel import OnClassModel
import obonet
from typing import Dict, Tuple
from tqdm import tqdm


## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "model": None,
    "reference_obs_targets": "cell_ontology_class",
    "input_layer": None,
    "reference_layer": None,
    "max_iter": 100,
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "cl_nlp_emb_file": "resources_test/annotation_test_data/ontology/cl.ontology.nlp.emb",
    "cl_ontology_file": "resources_test/annotation_test_data/ontology/cl.ontology",
    "cl_obo_file": "resources_test/annotation_test_data/ontology/cl.obo",
    "output_compression": "gzip"
}
meta = {"resources_dir": "src/annotate/onclass"}
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
def setup_logger():
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger

logger = setup_logger()

def map_celltype_to_ontology_id(cl_obo_file: str) -> Tuple[Dict[str, str], Dict[str, str]]:
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

def predict_input_data(model: OnClassModel,
                       input_matrix: np.array,
                       input_modality: ad.AnnData,
                       id_to_name: dict,
                       obs_prediction: str,
                       obs_probability: str) -> ad.AnnData:
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
    obs_prediction : str
        The obs key for the predicted cell type.
    obs_probability : str
        The obs key for the predicted cell type probability.
        
    Returns
    -------
    ad.AnnData
        The input data Anndata object with the predicted cell types saved in obs.
    """
    corr_test_feature = model.ProcessTestFeature(
        test_feature=input_matrix,
        test_genes=input_modality.var_names,
        log_transform=False,
    )
    onclass_pred = model.Predict(corr_test_feature, use_normalize=False, refine=True, unseen_ratio=-1.0)
    pred_label = [model.i2co[ind] for ind in onclass_pred[2]]
    pred_cell_type_label = [id_to_name[id] for id in pred_label]
    
    input_modality.obs[obs_prediction] = pred_cell_type_label
    input_modality.obs[obs_probability] = np.max(onclass_pred[1], axis=1) / onclass_pred[1].sum(1)
    return input_modality

def set_var_index(adata, var_name):
    adata.var.index = [re.sub("\\.[0-9]+$", "", s) for s in adata.var[var_name]]
    return adata

def main():
    
    if (not par["model"] and not par["reference"]) or (par["model"] and par["reference"]):
        raise ValueError("Make sure to provide either 'model' or 'reference', but not both.")
    
    logger.info("Reading input data")
    input_mudata = mu.read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]].copy()
    
    # Set var names to the desired gene name format (gene synbol, ensembl id, etc.)
    input_modality = set_var_index(input_modality, par["var_query_gene_names"]) if par["var_query_gene_names"] else input_modality
    input_matrix = input_modality.layers[par["input_layer"]].toarray() if par["input_layer"] else input_modality.X.toarray()

    id_to_name, name_to_id = map_celltype_to_ontology_id(par["cl_obo_file"])
    

    if par["model"]:
        logger.info("Predicting cell types using pre-trained model")
        model = OnClassModel(cell_type_nlp_emb_file=par["cl_nlp_emb_file"],
                             cell_type_network_file=par["cl_ontology_file"])
        
        model.BuildModel(use_pretrain=par["model"], ngene=None)
    
    
    elif par["reference"]:
        logger.info("Reading reference data")
        model = OnClassModel(cell_type_nlp_emb_file=par["cl_nlp_emb_file"],
                             cell_type_network_file=par["cl_ontology_file"])
        
        reference_mudata = mu.read_h5mu(par["reference"])
        reference_modality = reference_mudata.mod[par["modality"]].copy()

        reference_modality.var["gene_symbol"] = list(reference_modality.var.index)
        reference_modality.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference_modality.var["ensemblid"]]

        logger.info("Detecting common vars based on ensembl ids")
        common_ens_ids = list(set(reference_modality.var.index).intersection(set(input_modality.var.index)))

        logger.info("  reference n_vars: %i", reference_modality.n_vars)
        logger.info("  input n_vars: %i", input_modality.n_vars)
        logger.info("  intersect n_vars: %i", len(common_ens_ids))
        assert len(common_ens_ids) >= 100, "The intersection of genes is too small."

        reference_matrix = reference_modality.layers[par["reference_layer"]].toarray() if par["reference_layer"] else reference_modality.X.toarray()

        logger.info("Training a model from reference...")
        labels = reference_modality.obs[par["reference_obs_target"]].tolist()
        labels_cl = [name_to_id[label] for label in labels]
        _ = model.EmbedCellTypes(labels_cl)
        (
            corr_train_feature,
            _,
            corr_train_genes,
            _,
        ) = model.ProcessTrainFeature(
            train_feature=reference_matrix,
            train_label=labels_cl,
            train_genes=reference_modality.var_names,
            test_feature=input_matrix,
            test_genes=input_modality.var_names,
            log_transform=False,
        )
        model.BuildModel(ngene=len(corr_train_genes))
        model.Train(corr_train_feature,
                    labels_cl,
                    max_iter=par["max_iter"])
        
    
    logger.info(f"Predicting cell types")
    input_modality = predict_input_data(model,
                                        input_matrix,
                                        input_modality,
                                        id_to_name,
                                        par["output_obs_predictions"],
                                        par["output_obs_probability"])
    logger.info("Writing output data")
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])
    
if __name__ == "__main__":
    main()