import sys
from mudata import read_h5ad

## VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered.h5mu",
    "modality": "rna",
    "overwrite_existing_key": False,
    "input_obs_key": None,
    "output_obs_key": "index_copy",
    "output": "output.h5mu",
    "output_compression": "gzip",
}
meta = {"resources_dir": "src/metadata/copy_obs"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading from %s, modality %s", par["input"], par["modality"])
adata = read_h5ad(par["input"], mod=par["modality"])


def duplicate_obs(adata, input_key, output_key):
    if input_key:
        logger.info(f"Copying .obs key {input_key} to {output_key}")
        adata.obs[output_key] = adata.obs[input_key].copy()
    else:
        logger.info(f"Copying .obs index to {output_key}")
        adata.obs[output_key] = adata.obs.index.copy()


if par["output_obs_key"] not in adata.obs:
    duplicate_obs(adata, par["input_obs_key"], par["output_obs_key"])

else:
    if not par["overwrite_existing_key"]:
        raise ValueError(
            f"--output_obs_key already exists: `{par['output_obs_key']}`. Data can not be duplicated."
        )

    logger.warning(
        f"--output_obs_key already exists: `{par['output_obs_key']}`. Data in par['output_obs_key'] will be overwritten."
    )
    duplicate_obs(adata, par["input_obs_key"], par["output_obs_key"])


logger.info("Write output to mudata file")
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], adata, None
)
