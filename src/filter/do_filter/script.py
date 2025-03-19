import mudata as mu
import numpy as np
import sys

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "obs_filter": ["filter_none", "filter_with_random"],
    "var_filter": ["filter_with_random"],
    "output": "output.h5mu",
}

mdata = mu.read_h5mu(par["input"])
mdata.mod["rna"].obs["filter_none"] = np.repeat(True, mdata.mod["rna"].n_obs)
mdata.mod["rna"].obs["filter_with_random"] = np.random.choice(
    a=[False, True], size=mdata.mod["rna"].n_obs
)
mdata.mod["rna"].var["filter_with_random"] = np.random.choice(
    a=[False, True], size=mdata.mod["rna"].n_vars
)
mod = "rna"
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading %s, modality: %s", par["input"], par["modality"])
adata = mu.read_h5ad(par["input"], mod=par["modality"])

obs_filt = np.repeat(True, adata.n_obs)
var_filt = np.repeat(True, adata.n_vars)

par["obs_filter"] = par["obs_filter"] if par["obs_filter"] else []
par["var_filter"] = par["var_filter"] if par["var_filter"] else []

for obs_name in par["obs_filter"]:
    logger.info(
        "Filtering modality '%s' observations by .obs['%s']", par["modality"], obs_name
    )
    if obs_name not in adata.obs:
        raise ValueError(f".mod[{par['modality']}].obs[{obs_name}] does not exist.")
    if obs_name in adata.obs:
        obs_filt &= adata.obs[obs_name]

for var_name in par["var_filter"]:
    logger.info(
        "Filtering modality '%s' variables by .var['%s']", par["modality"], var_name
    )
    if var_name not in adata.var:
        raise ValueError(f".mod[{par['modality']}].var[{var_name}] does not exist.")
    if var_name in adata.var:
        var_filt &= adata.var[var_name]

adata_filtered = adata[obs_filt, var_filt]

logger.info("Writing h5mu to file %s.", par["output"])
write_h5ad_to_h5mu_with_compression(
    par["output"],
    par["input"],
    par["modality"],
    adata_filtered,
    par["output_compression"],
)
