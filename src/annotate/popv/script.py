import sys
import subprocess
import logging

import scanpy
import popv
from torch.cuda import is_available as cuda_is_available

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
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix_rna_anndata.h5ad",
    "output": "foo.h5ad",
    "compression": "gzip",
    # Component params
    "tissue": "Blood",
    "methods": ["scvi", "scanvi"],
    "obs_covariate": "batch",
    "obs_cell_type_key": "none",
    "obs_cell_type_unknown_label": "unknown",
    "reference_cell_type_key": "cell_ontology_class",
    "reference_covariate": ["donor", "method"]
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
    logger.info(result.stdout)
    logger.info(result.stderr)


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
    logger.info('Prepping input data for cell typing')
    
    # Query data
    adata_query = scanpy.read_h5ad(par["input"].strip())    
    adata_query.obs_names_make_unique()
    # Disclaimer: Hack around cross-species mouse-human, so be cautious
    adata_query.var_names = adata_query.var_names.str.upper()
    adata_query = check_nonnegative_integers(adata_query)

    # Reference data
    logger.info('Prepping Tabula Sapiens {} reference data for cell typing'.format(par["tissue"]))
    refdata_url = get_tabula_sapiens_reference_url(par["tissue"])
    local_refdata = 'tabula_sapiens_reference_{}.h5ad'.format(par["tissue"])
    download_tabula_sapiens_reference_h5ad(refdata_url, local_refdata)
    adata_reference = scanpy.read_h5ad(local_refdata)
    adata_reference = check_nonnegative_integers(adata_reference)
    ref_cell_types = adata_reference.obs['cell_type_tissue'].unique().to_list()
    logger.info('Reference data includes {} cell types {}'.format(len(ref_cell_types), ', '.join(ref_cell_types)))

    min_celltype_size = np.min(adata_reference.obs.groupby('cell_ontology_class').size())
    n_samples_per_label = np.max((min_celltype_size, 100))
    
    
    
    
if __name__ == "__main__":
    logger.info('PopV cell type annotation component')
    
    if not isinstance(par["methods"], list):
        par["methods"] = par["methods"].replace(" ", "").split(',')
    assert len(par["methods"]) >= 1, 'Please, specify at least one method for cell typing.'
    logger.info("Cell typing methods: {}".format(par["methods"]))

    if not par["obs_cell_type_key"] == 'none':
        assert isinstance(par["obs_cell_type_unknown_label"], str), 'Please, specify unknown cell type label in your .obs obs_cell_type_key.'
    else:
        par["obs_cell_type_key"] = None

    main()