import sys
import logging
import mudata as mu
import anndata as ad
import re
import numpy as np
import os
from OnClass.OnClassModel import OnClassModel
import obonet

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/tmp_TS_Blood_filtered.h5mu",
    "model": None,
    "reference_obs_targets": "cell_ontology_class",
    "input_layer": None,
    "reference_layer": None,
    "max_iter": 100,
    "output_obs_predictions": None,
    "output_obs_probability": None,
    "cl_nlp_emb_file": "cl.ontology.nlp.emb",
    "cl_ontology_file": "cl.ontology",
    "output_compression": "gzip"
}
meta = {}
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

def map_celltype_to_ontology_id(cl_obo_file):
    graph = obonet.read_obo(cl_obo_file)
    cl_id_to_name = {id_: data.get("name") for id_, data in graph.nodes(data=True)}
    return cl_id_to_name

def predict_input_data(model, input_matrix, input_modality, id_to_name, reference_obs_target, obs_prediction, obs_probability):
    corr_test_feature = model.ProcessTestFeature(
        test_feature=input_matrix,
        test_genes=input_modality.var_names,
        log_transform=False,
    )
    onclass_pred = model.Predict(corr_test_feature, use_normalize=False, refine=True, unseen_ratio=-1.0)
    pred_label = [model.i2co[ind] for ind in onclass_pred[2]]
    pred_cell_type_label = [id_to_name[id] for id in pred_label]
    
    input_modality.obs[f"{reference_obs_target}_{obs_prediction}"] = pred_cell_type_label
    input_modality.obs[f"{reference_obs_target}_{obs_probability}"] = np.max(onclass_pred[1], axis=1) / onclass_pred[1].sum(1)
    return input_modality

def main():
    input_mudata = mu.read_h5mu(par["input"])
    input_modality = input_mudata.mod[par["modality"]].copy()
    

    id_to_name = map_celltype_to_ontology_id(os.path.join(meta["resources_dir"], "cl.obo"))
    obs_predictions = par["output_obs_predictions"] if par["output_obs_predictions"] else ["pred"] * len(par["reference_obs_targets"])
    obs_probabilities = par["output_obs_probability"] if par["output_obs_probability"] else ["prob"] * len(par["reference_obs_targets"])
    
    if par["input_layer"]:
        input_matrix = input_modality.layers[par["input_layer"]]
    else:
        input_matrix = input_modality.X
        
    if par["reference"]:
        
        model = OnClassModel(cell_type_nlp_emb_file=par["cl_nlp_emb_file"],
                            cell_type_network_file=par["cl_ontology_file"])
        _ = model.EmbedCellTypes()
        
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

        if par["reference_layer"]:
            reference_matrix = reference_modality.layers[par["reference_layer"]]
        else:
            reference_matrix = reference_modality.X

        (
            corr_train_feature,
            _,
            corr_train_genes,
            _,
        ) = model.ProcessTrainFeature(
            train_feature=reference_matrix,
            train_label=reference_modality.obs[reference_obs_target],
            train_genes=reference_modality.var_names,
            test_feature=input_matrix,
            test_genes=input_modality.var_names,
            log_transform=False,
        )
            
        for reference_obs_target, obs_prediction, obs_probability in zip(par["reference_obs_targets"], obs_predictions, obs_probabilities):

            model.BuildModel(ngene=len(corr_train_genes))
            model.Train(corr_train_feature,
                        reference_modality.obs[reference_obs_target],
                        max_iter=par["max_iter"],
                        save_model=f"model_{reference_obs_target}")
            
            input_modality = predict_input_data(model,
                                                input_matrix,
                                                input_modality,
                                                id_to_name,
                                                reference_obs_target,
                                                obs_prediction,
                                                obs_probability)

    elif par["model"]:
        for model_path, reference_obs_target, obs_prediction, obs_probability in zip(par["model"], par["reference_obs_targets"], obs_predictions, obs_probabilities):
            model = OnClassModel(cell_type_nlp_emb_file=os.path.join(meta["resources_dir"], "cl.ontology.nlp.emb"),
                                cell_type_network_file=os.path.join(meta["resources_dir"], "cl.ontology"))
            model.BuildModel(use_pretrain=model_path)
        
            input_modality = predict_input_data(model,
                                                input_matrix,
                                                input_modality,
                                                id_to_name,
                                                reference_obs_target,
                                                obs_prediction,
                                                obs_probability)
    
    input_mudata.mod[par["modality"]] = input_modality
    input_mudata.write_h5mu(par["output"], compression=par["output_compression"])
    
if __name__ == "__main__":
    main()