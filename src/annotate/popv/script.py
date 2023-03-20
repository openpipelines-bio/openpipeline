import sys
import re
import logging
import numpy as np
import mudata as mu
import anndata as ad
import popv

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


## VIASH START
par = {
  # 'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
  'input': 'resources_test/concat_test_data/human_brain_3k_filtered_feature_bc_matrix_subset.h5mu',
  'modality': 'rna',
  'query_obs_covariates': ["batch"],
  'query_obs_cell_type_key': 'value',
  'query_obs_cell_type_unknown_label': 'unknown',
  'reference': "resources_test/annotation_test_data/tmp_TS_Blood_filtered.h5ad",
  'reference_obs_cell_type_key': "cell_ontology_class",
  'reference_obs_covariate': ["donor", "method"],
  'output': 'output.h5mu',
  'compression': 'gzip',
  'methods': ["knn_on_scvi", "scanvi"],
  "output_plots": "plots.pdf"
}
meta = {}
## VIASH END


use_gpu = cuda_is_available() or mps_is_available()
logger.info('GPU enabled? %s', use_gpu)

# changes:
# - input dataset is a mudata object
# - expect reference to be download apriori
# - CL ontology is now part of the docker container
# - switch reference to ensembl ids as soon as possible
# - find common variables by intersecting the ensembl genes
#   TODO: perform ortholog mapping if need be
# - TODO: turn static variables into args: https://github.com/czbiohub/PopV/blob/main/popv/preprocessing.py

# define static variables
input_obs_label = None
reference_obs_batch = "donor_assay"
reference_obs_label = "cell_ontology_class"

def main(par, meta):
    ### PREPROCESSING REFERENCE ###
    logger.info("### PREPROCESSING REFERENCE ###")
    
    # take a look at reference data
    logger.info("Reading reference data '%s'", par["reference"])
    reference = ad.read_h5ad(par["reference"])
    
    logger.info("Setting reference var index to Ensembl IDs")
    reference.var["gene_symbol"] = list(reference.var.index)
    reference.var.index = [re.sub("\\.[0-9]+$", "", s) for s in reference.var["ensemblid"]]

    logger.info("Detect number of samples per label")
    min_celltype_size = np.min(reference.obs.groupby(reference_obs_label).size())
    n_samples_per_label = np.max((min_celltype_size, 100))


    ### PREPROCESSING INPUT ###
    logger.info("### PREPROCESSING INPUT ###")
    logger.info("Reading '%s'", par["input"])
    input = mu.read_h5mu(par["input"])
    input_modality = input.mod[par['modality']]

    # subset with var column
    if par["input_var_subset"]:
        logger.info("Subset input with .var['%s']", par["input_var_subset"])
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
    adata = popv.preprocessing.Process_Query(
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
        save_path_trained_models="/pwd/output_popv", # TODO: fix
        n_samples_per_label=n_samples_per_label,
        # hardcoded values
        cl_obo_folder="/opt/popv_cl_ontology/",
        use_gpu=use_gpu
    ).adata
    
    # Lower mem footprint
    del input_modality_intersect
    del reference_intersect

    logger.info('Annotate data')
    popv.annotation.annotate_data(
        adata=adata.adata,
        methods=par["methods"],
        save_path="/pwd/output_annotate", # TODO: fix
    )
    
    # Lower mem footprint
    del adata_query_reference
    
    logger.info('Adding summary statistics...')
    predictions_fn = os.path.join(par["output_dir"], 'predictions.csv')
    predictions = pd.read_csv(predictions_fn, index_col = 0)
    for col in predictions.columns:
        adata_query.obs[col] = predictions.loc[adata_query.obs_names][col]
    if 'scvi' in par["methods"]:
        scvi_latent_space_fn = os.path.join(par["output_dir"], 'scvi_latent.csv')
        scvi_latent_space = pd.read_csv(scvi_latent_space_fn, index_col=0)
        adata_query.obsm['X_scvi'] = scvi_latent_space.loc[adata_query.obs_names]
    if 'scanvi' in par["methods"]:
        scanvi_latent_space_fn = os.path.join(par["output_dir"], 'scanvi_latent.csv')
        scanvi_latent_space = pd.read_csv(scanvi_latent_space_fn, index_col=0)
        adata_query.obsm['X_scanvi'] = scanvi_latent_space.loc[adata_query.obs_names]

    if 'scvi' in par["methods"]:
        logger.info('Re-calculating embedding...')
        sc.pp.neighbors(adata_query, use_rep="X_scvi")
        sc.tl.umap(adata_query, min_dist=0.3)

    logger.info('Storing cell annotated data...')
    adata_query.write(
        adata_query_cell_typed_file,
        compression=par["compression"]
        )
    
    if par["plots"]:
        logger.info('Creating agreement plots...')
        all_prediction_keys = [
        "popv_knn_on_bbknn_prediction",
        "popv_knn_on_scvi_online_prediction",
        "popv_knn_on_scvi_offline_prediction",
        "popv_scanvi_online_prediction",
        "popv_scanvi_offline_prediction",
        "popv_svm_prediction",
        "popv_rf_prediction",
        "popv_onclass_prediction",
        "popv_knn_on_scanorama_prediction"
        ]
        obs_keys = adata_query.obs.keys()
        pred_keys = [key for key in obs_keys if key in all_prediction_keys]
        popv.make_agreement_plots(
            adata_query,
            methods=pred_keys,
            popv_prediction_key = 'popv_prediction',
            save_folder=par["output_dir"]
            )
        
        logger.info('Creating frequency plots...')
        ax = adata_query.obs['popv_prediction_score'].value_counts().sort_index().plot.bar()
        ax.set_xlabel('Score')
        ax.set_ylabel("Frequency")
        ax.set_title("PopV Prediction Score Frequency")
        figpath = os.path.join(par["output_dir"], "prediction_score_barplot.pdf")
        ax.get_figure().savefig(figpath, bbox_inches="tight", dpi=300)
        
        ax = adata_query.obs.groupby('popv_prediction')['popv_prediction_score'].mean().plot.bar()
        ax.set_ylabel('Celltype')
        ax.set_xlabel('Averge number of methods in agreement')
        ax.set_title('Agreement per method by cell type')
        figpath = os.path.join(par["output_dir"], "percelltype_agreement_barplot.pdf")
        ax.get_figure().savefig(figpath, bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    logger.info('PopV cell type annotation component')
    
    input_file_name = ''.join(par["input"].split('/')[-1].split('.')[:-1])
    adata_query_cell_typed_file = '{}/{}_cell_typed.h5ad'.format(par["output_dir"], input_file_name)
    logger.info('PopV output can be found: {}'.format(par["output_dir"]))
    logger.info('Cell typed file stored: {}'.format(adata_query_cell_typed_file))

    if len(par["tissue_tabula_sapiens"]) >= 1:
        logger.info('Cell typing done on Tabula Sapiens reference data.')
        own_data = False
    elif len(par["tissue_reference_file"]) >= 1:
        logger.info('Cell typing done on user-provided provided reference data.')
        own_data = True
    else:
        raise BaseException("Please, tissue_tabula_sapiens or tissue_reference_file...")
    
    if not isinstance(par["methods"], list):
        par["methods"] = par["methods"].replace(" ", "").split(',')
    assert len(par["methods"]) >= 1, 'Please, specify at least one method for cell typing.'
    logger.info("Cell typing methods: {}".format(par["methods"]))

    if not par["query_obs_cell_type_key"] == 'none':
        assert isinstance(par["query_obs_cell_type_unknown_label"], str), 'Please, specify unknown cell type label in your .obs obs_cell_type_key.'
    else:
        par["query_obs_cell_type_key"] = None
        
    logger.info('Making file paths absolute')
    par["input"] = Path(par["input"].strip()).resolve()
    par["tissue_reference_file"] = Path(par["tissue_reference_file"].strip()).resolve()    

    main(par, meta)

