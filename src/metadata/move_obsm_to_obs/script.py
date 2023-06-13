from mudata import read_h5mu
from functools import partial
from pandas.errors import MergeError

import logging
from sys import stdout

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

### VIASH START
par = {
    "input": "work/f5/5f6365898ca5a42a360301a0c9e200/TSP15_Eye_ScleraEtc_10X_2_1.add_id.output.h5mu",
    "output": "foo.h5mu",
    "modality": "rna",
    "obsm_key": "sample_id",
}
### VIASH END



logger.info("Read mudata from file")
mdata = read_h5mu(par['input'])
try:
    mod_data = mdata.mod[par['modality']]
except KeyError:
    raise ValueError(f"Modality {par['modality']} does not exist.")

logger.info("Moving .obm key %s", par["obsm_key"])
try:
    obsm_matrix = mod_data.obsm[par["obsm_key"]].copy()
except KeyError:
    raise ValueError(f".obsm key {par['obsm_key']} was not found in "
                     f".obsm slot for modality {par['modality']}.")


obsm_matrix.rename(partial("{key}_{}".format, key=par["obsm_key"]), 
                   axis="columns", copy=False, inplace=True)

original_n_obs = len(mod_data.obs)
try:
    logger.info(f".obs names: {mod_data.obs_names}")
    logger.info(f".obsm index: {obsm_matrix.index}")
    mod_data.obs = mod_data.obs.merge(obsm_matrix, how="left",
                                      validate="one_to_one",
                                      left_index=True, right_index=True)
except MergeError as e:
    raise ValueError(f"Could not join .obsm matrix at {par['obsm_key']} to .obs because there "
                     "are some observation that are not overlapping between the two matrices "
                     "(indexes should overlap). This is either a bug or your mudata file is corrupt.")
del mod_data.obsm[par["obsm_key"]]

logger.info("Write output to mudata file")
mdata.write_h5mu(par['output'], compression=par["output_compression"])



