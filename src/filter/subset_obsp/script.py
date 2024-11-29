import sys
import mudata as mu

### VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu',
  'modality': 'rna',
  'input_obsp_key': 'distances',
  'input_obs_key': 'leiden',
  'input_obs_value': '1',
  'output_obsm_key': "leiden_1",
  'output': 'subset_obsp_output.h5mu',
  'output_compression': None,
}
### VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
logger = setup_logger()

def main():
    logger.info(f"Reading {par['input']}")
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]]

    logger.info(f"Subset columns of obsp matrix under {par['input_obsp_key']} based on {par['input_obs_key']} == {par['input_obs_value']}")
    # .obsp, .obs and .obsm index and .obsp columns all have a dimension length of `n_obs`
    # the index dimensions remain unaltered, but .obsp columns will be subset 
    obsp = adata.obsp[par["input_obsp_key"]]
    idx = adata.obs[par["input_obs_key"]].astype(str) == par["input_obs_value"]
    obsm_subset = obsp[:, idx]

    logger.info(f"Writing subset obsp matrix to .obsm {par['output_obsm_key']}")
    adata.obsm[par["output_obsm_key"]] = obsm_subset

    logger.info(f"Writing output to {par['output']}")
    mdata.write_h5mu(par["output"], compression=par["output_compression"])


if __name__ == '__main__':
    main()
