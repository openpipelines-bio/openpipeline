import sys
import pandas as pd
from mudata import read_h5ad

### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "uns_key": "metrics_cellranger",
    "output": "foo.h5mu",
    "modality": "rna",
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Read modality %s from file %s", par["modality"], par["input"])
mod_data = read_h5ad(par["input"], mod=par["modality"])

logger.info("Joining uns to obs")
# get data frame
uns_df = mod_data.uns[par["uns_key"]]

# check for overlapping colnames
intersect_keys = uns_df.keys().intersection(mod_data.obs.keys())
obs_drop = mod_data.obs.drop(intersect_keys, axis=1)

# create data frame to join
uns_df_rep = uns_df.loc[uns_df.index.repeat(mod_data.n_obs)]
uns_df_rep.index = mod_data.obs_names

# create new obs
mod_data.obs = pd.concat([obs_drop, uns_df_rep], axis=1)

logger.info("Write output to mudata file")
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], mod_data, par["output_compression"]
)
