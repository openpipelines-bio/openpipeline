import sys
import re
import tempfile
import typing
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

# where to find the obo files
cl_obo_folder = "/opt/PopV/resources/ontology/"

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    # "input": "resources_test/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu",
    "modality": "rna",
    "reference": "resources_test/annotation_test_data/tmp_TS_Blood_filtered.h5ad",
    "input_obs_batch": None,
    "input_layer": None,
    "input_obs_label": None,
    "input_var_subset": None,
    "unknown_celltype_label": "unknown",
    "reference_layer": None,
    "reference_obs_label": "cell_ontology_class",
    "reference_obs_batch": "donor_assay",
    "output": "output.h5mu",
    "output_compression": "gzip",
    "methods": [
        # "celltypist",
        # "knn_on_bbknn",
        # "knn_on_scanorama",
        # "knn_on_scvi",
        "rf",
        # "scanvi",
        "svm",
    ]
}
meta = {}
# for debugging the obo folder can be somewhere local
cl_obo_folder = "popv_cl_ontology/"
## VIASH END

sys.path.append(meta["resources_dir"])
# START TEMPORARY WORKAROUND setup_logger
# reason: resources aren't available when using Nextflow fusion
# from setup_logger import setup_logger
def setup_logger():
    import logging
    from sys import stdout

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    console_handler = logging.StreamHandler(stdout)
    logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
    console_handler.setFormatter(logFormatter)
    logger.addHandler(console_handler)

    return logger
# END TEMPORARY WORKAROUND setup_logger
logger = setup_logger()

use_gpu = cuda_is_available()
logger.info("GPU enabled? %s", use_gpu)

# Helper functions
def get_X(adata: ad.AnnData, layer: typing.Optional[str], var_index: typing.Optional[str]):
    """Fetch the counts data from X or a layer. Subset columns by var_index if so desired."""
    if var_index:
        adata = adata[:, var_index]
    if layer:
        return adata.layers[layer]
    else:
        return adata.X
def get_obs(adata: ad.AnnData, obs_par_names):
    """Subset the obs dataframe to just the columns defined by the obs_label and obs_batch."""
    obs_columns = [par[x] for x in obs_par_names if par[x]]
    return adata.obs[obs_columns]
def get_var(adata: ad.AnnData, var_index: list[str]):
    """Fetch the var dataframe. Subset rows by var_index if so desired."""
    return adata.var.loc[var_index]

def main(par, meta):
    assert len(par["methods"]) >= 1, "Please, specify at least one method for cell typing."
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
    input_modality = input.mod[par["modality"]]

    # subset with var column
    if par["input_var_subset"]:
        logger.info("Subset input with .var['%s']", par["input_var_subset"])
        assert par["input_var_subset"] in input_modality.var, f"--input_var_subset='{par['input_var_subset']}' needs to be a column in .var"
        input_modality = input_modality[:,input_modality.var[par["input_var_subset"]]]

    ### ALIGN REFERENCE AND INPUT ###
    logger.info("### ALIGN REFERENCE AND INPUT ###")

    logger.info("Detecting common vars based on ensembl ids")
    common_ens_ids = list(set(reference.var.index).intersection(set(input_modality.var.index)))
    
    logger.info("  reference n_vars: %i", reference.n_vars)
    logger.info("  input n_vars: %i", input_modality.n_vars)
    logger.info("  intersect n_vars: %i", len(common_ens_ids))
    assert len(common_ens_ids) >= 100, "The intersection of genes is too small."

    # subset input objects to make sure popv is using the data we expect
    input_modality = ad.AnnData(
        X = get_X(input_modality, par["input_layer"], common_ens_ids),
        obs = get_obs(input_modality, ["input_obs_label", "input_obs_batch"]),
        var = get_var(input_modality, common_ens_ids)
    )
    reference = ad.AnnData(
        X = get_X(reference, par["reference_layer"], common_ens_ids),
        obs = get_obs(reference, ["reference_obs_label", "reference_obs_batch"]),
        var = get_var(reference, common_ens_ids)
    )

    # remove layers that 
    
    ### ALIGN REFERENCE AND INPUT ###
    logger.info("### ALIGN REFERENCE AND INPUT ###")

    with tempfile.TemporaryDirectory(prefix="popv-", dir=meta["temp_dir"]) as temp_dir:
        logger.info("Run PopV processing")
        pq = popv.preprocessing.Process_Query(
            # input
            query_adata=input_modality,
            query_labels_key=par["input_obs_label"],
            query_batch_key=par["input_obs_batch"],
            query_layers_key=None, # this is taken care of by subset
            # reference
            ref_adata=reference,
            ref_labels_key=par["reference_obs_label"],
            ref_batch_key=par["reference_obs_batch"],
            # options
            unknown_celltype_label=par["unknown_celltype_label"],
            n_samples_per_label=n_samples_per_label,
            # pretrained model
            # Might need to be parameterized at some point
            prediction_mode="retrain",
            pretrained_scvi_path=None,
            # outputs
            # Might need to be parameterized at some point
            save_path_trained_models=temp_dir,
            # hardcoded values
            cl_obo_folder=cl_obo_folder,
            accelerator='cuda' if use_gpu else None
        )
        method_kwargs = {}
        if 'scanorama' in par['methods']:
            method_kwargs['scanorama'] = {'approx': False}
        logger.info("Annotate data")
        popv.annotation.annotate_data(
            adata=pq.adata,
            methods=par["methods"],
            methods_kwargs=method_kwargs
        )

    popv_input = pq.adata[input_modality.obs_names]

    # select columns starting with "popv_"
    popv_obs_cols = popv_input.obs.columns[popv_input.obs.columns.str.startswith("popv_")]

    # create new data frame with selected columns
    df_popv = popv_input.obs[popv_obs_cols]

    # remove prefix from column names
    df_popv.columns = df_popv.columns.str.replace("popv_", "")

    # store output in mudata .obsm
    input.mod[par["modality"]].obsm["popv_output"] = df_popv

    # copy important output in mudata .obs
    for col in ["popv_prediction"]:
        if col in popv_input.obs.columns:
            input.mod[par["modality"]].obs[col] = popv_input.obs[col]

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

