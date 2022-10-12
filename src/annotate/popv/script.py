import os
import sys
import subprocess
import logging

import pandas as pd
import matplotlib
import scanpy
import popv
from torch.cuda import is_available as cuda_is_available
import numpy as np

try:
    from torch.backends.mps import is_available as mps_is_available
except ModuleNotFoundError:
    # Older pytorch versions
    # MacOS GPUs
    def mps_is_available():
        return False


logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(sys.stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)


### VIASH START
par = {
    # I/O params
    "input": "./resources_test/annotation_test_data/blood_test_query.h5ad",
    "output_dir": "./resources_test/annotation_test_data/",
    "compression": "gzip",
    # Component params
    "tissue_tabula_sapiens": "Blood",
    "tissue_reference_file": "./resources_test/annotation_test_data/blood_test_reference.h5ad",
    "methods": ["scvi", "scanvi"],
    "query_obs_covariate": ["batch"],
    "query_obs_cell_type_key": "none",
    "query_obs_cell_type_unknown_label": "unknown",
    "reference_obs_cell_type_key": "cell_ontology_class",
    "reference_obs_covariate": ["donor", "method"],
    "ontology_files_path": "./resources_test/annotation_test_data/popv_cl_ontology/",
    "plots": False
    }
### VIASH END


def get_tabula_sapiens_reference_url(tissue):
    if tissue == 'Bladder':
        refdata_url ='https://www.dropbox.com/s/p5x1lb0jyl8293c/Bladder.h5ad?dl=1'
    elif tissue == 'Blood':
        refdata_url = 'https://www.dropbox.com/s/4cg6zj340oelhlg/Blood.h5ad?dl=1'
    elif tissue == 'Bone_Marrow':
        refdata_url = 'https://www.dropbox.com/s/rwfovoyafpd64io/Bone_Marrow.h5ad?dl=1'
    elif tissue == 'Fat':
        refdata_url = 'https://www.dropbox.com/s/if1d7iloovi8e9o/Fat.h5ad?dl=1'
    elif tissue == 'Heart':
        refdata_url = 'https://www.dropbox.com/s/0udrdzjl2z087jj/Heart.h5ad?dl=1'
    elif tissue == 'Kidney':
        refdata_url = 'https://www.dropbox.com/s/8sx9fhjfgnyjgdz/Kidney.h5ad?dl=1'
    elif tissue == 'Large_Intestine':
        refdata_url = 'https://www.dropbox.com/s/272sajn0hkj62le/Large_Intestine.h5ad?dl=1'
    elif tissue == 'Liver':
        refdata_url = 'https://www.dropbox.com/s/g0ahumalnm0mp38/Liver.h5ad?dl=1'
    elif tissue == 'Lung':
        refdata_url = 'https://www.dropbox.com/s/2kuzdamjevev2ci/Lung.h5ad?dl=1'
    elif tissue == 'Lymph_Node':
        refdata_url = 'https://www.dropbox.com/s/tetuh62010uothb/Lymph_Node.h5ad?dl=1'
    elif tissue == 'Mammary':
        refdata_url = 'https://www.dropbox.com/s/krm4pv4ev6cynns/Mammary.h5ad?dl=1'
    elif tissue == 'Muscle':
        refdata_url = 'https://www.dropbox.com/s/0jhvnoy49rvrlqn/Muscle.h5ad?dl=1'
    elif tissue == 'Pancreas':
        refdata_url = 'https://www.dropbox.com/s/kn0zodnmxwx0yhe/Pancreas.h5ad?dl=1'
    elif tissue == 'Prostate':
        refdata_url = 'https://www.dropbox.com/s/040fb5jr0zcur7h/Prostate.h5ad?dl=1'
    elif tissue == 'Salivary Gland':
        refdata_url = 'https://www.dropbox.com/s/rwia1ji7eztga6b/Salivary_Gland.h5ad?dl=1'
    elif tissue == 'Skin':  
        refdata_url = 'https://www.dropbox.com/s/ucvdksq2jnug2nh/Skin.h5ad?dl=1'
    elif tissue == 'Small_Intestine':
        refdata_url = 'https://www.dropbox.com/s/06ia5n2yex3dq8j/Small_Intestine.h5ad?dl=1'
    elif tissue == 'Spleen':
        refdata_url = 'https://www.dropbox.com/s/m2d0gme847qdhr1/Spleen.h5ad?dl=1'
    elif tissue == 'Thymus':
        refdata_url = 'https://www.dropbox.com/s/i84bcyk87scesml/Thymus.h5ad?dl=1'
    elif tissue == 'Trachea':
        refdata_url = 'https://www.dropbox.com/s/ppt7b6w73gvceap/Trachea.h5ad?dl=1'
    elif tissue == 'Vasculature':
        refdata_url = 'https://www.dropbox.com/s/1eq0zamel5etmoq/Vasculature.h5ad?dl=1'
    return refdata_url


def download_tabula_sapiens_reference_h5ad(refdata_url, local_refdata):
    logger.info('Downloading reference data')
    result = subprocess.run(
        ["wget", refdata_url, "-O", local_refdata], capture_output=True, text=True
        )
    # logger.info(result.stdout)
    # logger.info(result.stderr)


def check_nonnegative_integers(adata):
    if popv._check_nonnegative_integers(adata.X):
        logger.info('Raw counts present in .X')
    elif popv._check_nonnegative_integers(adata.raw.X):
        logger.info('Raw counts present in .raw.X')
        adata.X = adata.raw.X
    else:
        raise BaseException('Make sure data .X or .raw.X contains raw counts.')
    return adata


def main():
    
    logger.info('Prepping query data for cell typing')
    adata_query = scanpy.read_h5ad(par["input"].strip())
    # TODO: Hack around cross-species mouse-human, so be cautious
    # TODO: implement cross-species mapping
    adata_query.var_names = adata_query.var_names.str.upper()
    adata_query.obs_names_make_unique()
    adata_query.var_names_make_unique()
    adata_query = check_nonnegative_integers(adata_query)

    # TODO: add possibility to come up with own reference data e.g. atlasses
    if own_data:
        logger.info('Prepping user-provided {} reference data for cell typing'.format(par["tissue_reference_file"].strip()))
        adata_reference = scanpy.read_h5ad(par["tissue_reference_file"].strip())
    else:
        logger.info('Prepping Tabula Sapiens {} reference data for cell typing'.format(par["tissue_tabula_sapiens"]))
        refdata_url = get_tabula_sapiens_reference_url(par["tissue_tabula_sapiens"])
        local_refdata = 'tabula_sapiens_reference_{}.h5ad'.format(par["tissue_tabula_sapiens"])
        download_tabula_sapiens_reference_h5ad(refdata_url, local_refdata)
        adata_reference = scanpy.read_h5ad(local_refdata)
    adata_reference.obs_names_make_unique()
    adata_reference.var_names_make_unique()
    adata_reference = check_nonnegative_integers(adata_reference)
    ref_cell_types = adata_reference.obs['cell_type_tissue'].unique().to_list()
    logger.info('Reference data includes {} cell types {}'.format(len(ref_cell_types), ', '.join(ref_cell_types)))

    min_celltype_size = np.min(adata_reference.obs.groupby('cell_ontology_class').size())
    n_samples_per_label = np.max((min_celltype_size, 100))
    
    common_markers = list(set(adata_reference.var.index).intersection(set(adata_query.var.index)))
    logger.info(('Common features between query and reference data: {}'.format(len(common_markers))))
    adata_query_intersection = adata_query[:, common_markers]
    adata_reference_intersection = adata_reference[:, common_markers]
    assert adata_query_intersection.var.index.shape == adata_reference_intersection.var.index.shape, 'Something went wrong in the process of finding common markers!'
    
    # Lower mem footprint
    del adata_reference
    
    # TODO: add highly variable genes param to component
    # TODO: add training_mode to component
    # TODO: add pretrained_scvi option
    # TODO: add pretrained_scANVI option
    logger.info('PopV process query...')
    adata_query_reference = popv.process_query(
        query_adata=adata_query_intersection,
        ref_adata=adata_reference_intersection,
        save_folder=par["output_dir"],
        query_batch_key='batch',
        query_labels_key=None,
        unknown_celltype_label='unknown',
        pretrained_scvi_path=None,
        training_mode='offline',
        ref_labels_key="cell_ontology_class",
        ref_batch_key=['donor', 'method'],
        n_samples_per_label=n_samples_per_label,
        hvg=True)
    
    # Lower mem footprint
    del adata_query_intersection
    del adata_reference_intersection

    # TODO: assess CPU performance
    # TODO: assess GPU performance
    # TODO: add pretrained_scvi option
    # TODO: add pretrained_scANVI option
    logger.info('PopV annotate data...')
    popv.annotate_data(
        adata=adata_query_reference,
        methods=par["methods"],
        save_path=par["output_dir"],
        pretrained_scvi_path=None,
        pretrained_scanvi_path=None,
        onclass_ontology_file="{}/cl.ontology".format(par["ontology_files_path"]),
        onclass_obo_fp="{}/cl.obo".format(par["ontology_files_path"]),
        onclass_emb_fp="{}/cl.ontology.nlp.emb".format(par["ontology_files_path"]),
        scvi_max_epochs=None,
        scanvi_max_epochs=None,
        use_gpu=(cuda_is_available() or mps_is_available())
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
    adata_query.write(adata_query_cell_typed_file, compression='gzip')
    
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
    
    logger.info('GPU enabled? {}'.format((cuda_is_available() or mps_is_available())))
    
    input_file_name = par["input"].split('/')[-1].split('.')[-1]
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

    main()