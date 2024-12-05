import mudata as mu
import numpy as np
import sys
from operator import le, ge
from pandas.api.types import is_float_dtype


### VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "var_name_filter": "filter_with_counts",
    "min_fraction": 0,
    "max_fraction": 1,
    "output_compression": "gzip",
}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()

logger.info("Reading input data")
mdata = mu.read_h5mu(par["input"])

mdata.var_names_make_unique()

mod = par["modality"]
logger.info("Processing modality %s.", mod)
data = mdata.mod[mod]

logger.info("\tUnfiltered data: %s", data)

logger.info("\tComputing aggregations.")


def apply_filter_to_mask(mask, base, filter, comparator):
    new_filt = np.ravel(comparator(base, filter))
    num_removed = np.sum(np.invert(new_filt) & mask)
    mask &= new_filt
    return num_removed, mask


try:
    fraction = data.obs[par["obs_fraction_column"]]
except KeyError:
    raise ValueError(f"Could not find column '{par['obs_fraction_column']}'")
if not is_float_dtype(fraction):
    raise ValueError(
        f"Column '{par['obs_fraction_column']}' does not contain float datatype."
    )
if fraction.max() > 1:
    raise ValueError(f"Column '{par['obs_fraction_column']}' contains values > 1.")
if fraction.min() < 0:
    raise ValueError(f"Column '{par['obs_fraction_column']}' contains values < 0.")


# Filter cells
filters = (
    (
        "min_fraction",
        fraction,
        ge,
        "\tRemoving %s cells with <%s percentage mitochondrial reads.",
    ),
    (
        "max_fraction",
        fraction,
        le,
        "\tRemoving %s cells with >%s percentage mitochondrial reads.",
    ),
)

keep_cells = np.repeat(True, data.n_obs)
for filter_name_or_value, base, comparator, message in filters:
    try:
        filter = par[filter_name_or_value]
    except KeyError:
        filter = filter_name_or_value
    if filter is not None:
        num_removed, keep_cells = apply_filter_to_mask(
            keep_cells, base, filter, comparator
        )
        logger.info(message, num_removed, filter)

data.obs[par["obs_name_filter"]] = keep_cells

logger.info("\tFiltered data: %s", data)
logger.info("Writing output data to %s", par["output"])
mdata.write_h5mu(par["output"], compression=par["output_compression"])

logger.info("Finished")
