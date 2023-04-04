
import sys
import re
import logging
import numpy as np
import mudata as mu
import anndata as ad
import popv

# todo: is this still needed?
from torch.cuda import is_available as cuda_is_available
try:
    from torch.backends.mps import is_available as mps_is_available
except ModuleNotFoundError:
    # Older pytorch versions
    # MacOS GPUs
    def mps_is_available():
        return False

# set up logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

# where to find the obo files
cl_obo_folder = "/opt/popv_cl_ontology/"

## VIASH START
par = {
    'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
    # 'input': 'resources_test/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu',
    'modality': 'rna',
    'reference': 'resources_test/annotation_test_data/tmp_TS_Blood_filtered.h5ad',
    'input_obs_batch': None,
    'input_obs_layer': None,
    'input_obs_labels': None,
    'input_var_subset': None,
    'unknown_celltype_label': 'unknown',
    'reference_obs_labels': 'cell_ontology_class',
    'reference_obs_batch': 'donor_assay',
    'reference_models': None,
    'output': 'output.h5mu',
    'output_compression': 'gzip',
    'output_models': 'output.tar.gz',
    'methods': [
        # 'celltypist',
        # 'knn_on_bbknn',
        # 'knn_on_scanorama',
        # 'knn_on_scvi',
        'rf',
        # 'scanvi',
        'svm',
    ],
    'prediction_mode': None
}
meta = {}
# for debugging the obo folder can be somewhere local
cl_obo_folder = "popv_cl_ontology/"
## VIASH END

use_gpu = cuda_is_available() or mps_is_available()
logger.info('GPU enabled? %s', use_gpu)

# changes:
# - input dataset is a mudata object
# - expect reference to be download apriori
# - CL ontology is now part of the docker container
# - switch reference to ensembl ids as soon as possible
# - find common variables by intersecting the ensembl genes
# - store output in obsm, and store some values in the obs
#   TODO: perform ortholog mapping if need be

def main(par, meta):
    assert len(par["methods"]) >= 1, 'Please, specify at least one method for cell typing.'
    logger.info("Cell typing methods: {}".format(par["methods"]))

    ### PREPROCESSING REFERENCE ###
    logger.info("### PREPROCESSING REFERENCE ###")
    
    # take a look at reference data
    logger.info("Reading reference data '%s'", par["reference"])
    reference = ad.read_h5ad(par["reference"])
    
    logger.info("Setting reference var index to Ensembl IDs")
    reference.var["gene_symbol"] = list(reference.var.index)
    reference.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference.var["ensemblid"]]

    logger.info("Detect number of samples per label")
    min_celltype_size = np.min(reference.obs.groupby(par["reference_obs_batch"]).size())
    n_samples_per_label = np.max((min_celltype_size, 100))

    ### PREPROCESSING INPUT ###
    logger.info("### PREPROCESSING INPUT ###")
    logger.info("Reading '%s'", par["input"])
    input = mu.read_h5mu(par["input"])
    input_modality = input.mod[par['modality']]

    # subset with var column
    if par["input_var_subset"]:
        logger.info("Subset input with .var['%s']", par["input_var_subset"])
        assert par["input_var_subset"] in input_modality.var, f"--input_var_subset='{par['input_var_subset']}' needs to be a column in .var"
        input_modality = input_modality[:,input_modality.var[par["input_var_subset"]]].copy()

    ### ALIGN REFERENCE AND INPUT ###
    logger.info("### ALIGN REFERENCE AND INPUT ###")

    logger.info("Detecting common vars based on ensembl ids")
    common_ens_ids = list(set(reference.var.index).intersection(set(input_modality.var.index)))
    
    logger.info('  reference n_vars: %i', reference.n_vars)
    logger.info('  input n_vars: %i', input_modality.n_vars)
    logger.info('  intersect n_vars: %i', len(common_ens_ids))
    assert len(common_ens_ids) >= 100, "The intersection of genes is too small."

    # perform intersection
    input_modality = input_modality[:, common_ens_ids].copy()
    reference = reference[:, common_ens_ids].copy()
    
    ### ALIGN REFERENCE AND INPUT ###
    logger.info("### ALIGN REFERENCE AND INPUT ###")
    # TODO: add training_mode to component
    logger.info('Run PopV processing')
    pq = popv.preprocessing.Process_Query(
        # input
        query_adata=input_modality,
        query_labels_key=par["input_obs_labels"],
        query_batch_key=par["input_obs_batch"],
        query_layers_key=par["input_obs_layer"],
        # reference
        ref_adata=reference,
        ref_labels_key=par["reference_obs_labels"],
        ref_batch_key=par["reference_obs_batch"],
        # pretrained model
        prediction_mode='retrain',
        pretrained_scvi_path=None, # TODO: implement
        # options
        unknown_celltype_label='unknown',
        # outputs
        save_path_trained_models="output_popv", # TODO: fix output path
        n_samples_per_label=n_samples_per_label,
        # hardcoded values
        cl_obo_folder=cl_obo_folder,
        use_gpu=use_gpu
    )
    # pq_adata_orig = pq.adata.copy()

    logger.info('Annotate data')
    popv.annotation.annotate_data(
        adata=pq.adata,
        methods=par["methods"]
    )

    popv_input = pq.adata[input.obs_names]

    # select columns starting with "popv_"
    popv_obs_cols = popv_input.obs.columns[popv_input.obs.columns.str.startswith('popv_')]

    # create new data frame with selected columns
    df_popv = popv_input.obs[popv_obs_cols]

    # remove prefix from column names
    df_popv.columns = df_popv.columns.str.replace('popv_', '')

    # store output in mudata .obsm
    input.mod[par['modality']].obsm["popv_output"] = df_popv

    # copy important output in mudata .obs
    for col in ["popv_prediction"]:
        if col in popv_input.obs.columns:
            input.mod[par['modality']].obs[col] = popv_input.obs[col]

    # code to explore how the output differs from the original
    # for attr in ["obs", "var", "uns", "obsm", "layers", "obsp"]:
    #     old_keys = set(getattr(pq_adata_orig, attr).keys())
    #     new_keys = set(getattr(pq.adata, attr).keys())
    #     diff_keys = list(new_keys.difference(old_keys))
    #     diff_keys.sort()
    #     print(f"{attr}:", flush=True)
    #     for key in diff_keys:
    #         print(f"  {key}", flush=True)
    
    # write output
    logger.info("Writing %s", par["output"])
    input.write_h5mu(par["output"], compression=par["output_compression"])

if __name__ == "__main__":
    main(par, meta)

