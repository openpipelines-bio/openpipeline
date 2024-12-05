import sys
from functools import partial
from pandas.errors import MergeError
from mudata import read_h5mu

## VIASH START
par = {
    "input": "input.h5mu",
    "modality": "mod1",
    "obsm_key": "obsm_key",
    "output": "output.h5mu",
    "output_compression": "gzip",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info("Read mudata from file")
mdata = read_h5mu(par["input"])
try:
    mod_data = mdata.mod[par["modality"]]
except KeyError:
    raise ValueError(f"Modality {par['modality']} does not exist.")

logger.info("Moving .obm key %s", par["obsm_key"])
try:
    obsm_matrix = mod_data.obsm[par["obsm_key"]].copy()
except KeyError:
    raise ValueError(
        f".obsm key {par['obsm_key']} was not found in "
        f".obsm slot for modality {par['modality']}."
    )


obsm_matrix.rename(
    partial("{key}_{}".format, key=par["obsm_key"]),
    axis="columns",
    copy=False,
    inplace=True,
)

original_n_obs = len(mod_data.obs)
try:
    logger.info(f".obs names: {mod_data.obs_names}")
    logger.info(f".obsm index: {obsm_matrix.index}")
    new_obs = mod_data.obs.drop(obsm_matrix.columns, axis=1, errors="ignore")
    new_obs = new_obs.merge(
        obsm_matrix,
        how="left",
        validate="one_to_one",
        left_index=True,
        right_index=True,
    )
    mod_data.obs = new_obs
except MergeError:
    raise ValueError(
        f"Could not join .obsm matrix at {par['obsm_key']} to .obs because there "
        "are some observation that are not overlapping between the two matrices "
        "(indexes should overlap). This is either a bug or your mudata file is corrupt."
    )
del mod_data.obsm[par["obsm_key"]]

logger.info("Write output to mudata file")
mdata.write_h5mu(par["output"], compression=par["output_compression"])
