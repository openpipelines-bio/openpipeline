import sys
import mudata as mu

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "modality": "rna",
    "input_obsp_key": "distances",
    "input_obs_key": "leiden",
    "input_obs_value": "1",
    "output_obsm_key": "leiden_1",
    "output": "subset_obsp_output.h5mu",
    "output_compression": None,
}
### VIASH END
sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def main():
    logger.info("Reading %s, modality", par["input"], par["modality"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])

    logger.info(
        f"Subset columns of obsp matrix under {par['input_obsp_key']} based on {par['input_obs_key']} == {par['input_obs_value']}"
    )
    # .obsp, .obs and .obsm index and .obsp columns all have a dimension length of `n_obs`
    # the index dimensions remain unaltered, but .obsp columns will be subset
    obsp = adata.obsp[par["input_obsp_key"]]
    idx = adata.obs[par["input_obs_key"]].astype(str) == par["input_obs_value"]
    # A Series object cannot be used as an indexer for a scipy sparse array
    # when the data type is a pandas boolean extension array because
    # extension arrays do not define .nonzero()
    # See https://github.com/pandas-dev/pandas/issues/46025
    idx = idx.to_numpy(dtype="bool", na_value=False)
    obsm_subset = obsp[:, idx]

    logger.info(f"Writing subset obsp matrix to .obsm {par['output_obsm_key']}")
    adata.obsm[par["output_obsm_key"]] = obsm_subset

    logger.info(
        "Writing output to %s, modality %s", par["output"], par["output_compression"]
    )
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], adata, par["output_compression"]
    )


if __name__ == "__main__":
    main()
