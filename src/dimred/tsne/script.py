import scanpy as sc
import mudata as mu
import sys
import anndata as ad

## VIASH START
par = {
    'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu',
    'modality': 'rna',
    'use_rep': 'X_pca',
    'output': 'output.h5mu',
    'output_compression': 'gzip',
    'obsm_output': 'X_tsne',
    'n_pcs': 50,
    'perplexity': 30,
    'min_dist': 0.5,
    'metric': 'euclidean',
    'early_exaggeration': 12,
    'learning_rate': 1000,
    'random_state': 0,
    'n_jobs': None,
}
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

logger.info("Reading %s", par["input"])
mdata = mu.read_h5mu(par["input"])

logger.info("Computing tSNE for modality '%s'", par['modality'])
data = mdata.mod[par['modality']]

if par['use_rep'] not in data.obsm.keys():
    raise ValueError(f"'{par['use_rep']}' was not found in .mod['{par['modality']}'].obsm. No precomputed PCA provided. Please run PCA first.")
temp_obsm = {par["use_rep"]: data.obsm[par["use_rep"]]}

# create temporary AnnData
# ... because sc.tl.umap doesn't allow to choose
# the obsm output slot
# ... also we can see scanpy is a data format dependency hell
temp_adata = ad.AnnData(
    obsm=temp_obsm,
    shape=data.shape
)

sc.tl.tsne(
    adata=temp_adata,
    n_pcs=par["n_pcs"],
    use_rep=par["use_rep"],
    perplexity=par["perplexity"],
    metric=par["metric"],
    early_exaggeration=par["early_exaggeration"],
    learning_rate=par["learning_rate"],
    random_state=par["random_state"],
    n_jobs=par["n_jobs"]
)

logger.info(f"Writing tSNE embeddings to .mod[{par['modality']}].obsm[{par['obsm_output']}]")
data.obsm[par['obsm_output']] = temp_adata.obsm['X_tsne']

logger.info(f"Writing tSNE metadata to .mod[{par['modality']}].uns['tsne']")
data.uns['tsne'] = temp_adata.uns['tsne']

logger.info("Writing to %s.", par["output"])
mdata.write_h5mu(filename=par["output"], compression=par["output_compression"])

logger.info("Finished")