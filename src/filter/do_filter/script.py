import mudata as mu
import numpy as np
import logging
from sys import stdout

### VIASH START
par = {
  'input': 'resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu',
  'modality': 'rna',
  'obs_filter': ['filter_none', 'filter_with_random'],
  'var_filter': ['filter_with_random'],
  'output': 'output.h5mu'
}

mdata = mu.read_h5mu(par["input"])
mdata.mod['rna'].obs["filter_none"] = np.repeat(True, mdata.mod['rna'].n_obs)
mdata.mod['rna'].obs["filter_with_random"] = np.random.choice(a=[False, True], size=mdata.mod['rna'].n_obs)
mdata.mod['rna'].var["filter_with_random"] = np.random.choice(a=[False, True], size=mdata.mod['rna'].n_vars)
mod = 'rna'
### VIASH END

logger = logging.getLogger()
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler(stdout)
logFormatter = logging.Formatter("%(asctime)s %(levelname)-8s %(message)s")
console_handler.setFormatter(logFormatter)
logger.addHandler(console_handler)

logger.info("Reading %s", par['input'])
mdata = mu.read_h5mu(par["input"])

mod = par["modality"]
logger.info("Processing modality '%s'", mod)

obs_filt = np.repeat(True, mdata.mod[mod].n_obs)
var_filt = np.repeat(True, mdata.mod[mod].n_vars)

for obs_name in par["obs_filter"]:
    logger.info("Filtering modality '%s' observations by .obs['%s']", mod, obs_name)
    if obs_name in mdata.mod[mod].obs:
        obs_filt &= mdata.mod[mod].obs[obs_name]
    else:
        logger.warning(".mod['%s'].obs['%s'] does not exist. Skipping.", mod, obs_name)

for var_name in par["var_filter"]:
    logger.info("Filtering modality '%s' variables by .var['%s']", mod, obs_name)
    if var_name in mdata.mod[mod].var:
        var_filt &= mdata.mod[mod].var[var_name]
    else:
        logger.warning(".mod['%s'.var['%s'] does not exist. Skipping.", mod, obs_name)

mdata.mod[mod] = mdata.mod[mod][obs_filt, var_filt].copy()

logger.info("Writing h5mu to file %s.", par["output"])
mdata.write_h5mu(par["output"], compression="gzip")
