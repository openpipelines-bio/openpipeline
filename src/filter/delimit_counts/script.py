import mudata as mu
import numpy as np
import sys
from operator import le, ge
from pandas.api.types import is_numeric_dtype


### VIASH START
par = {
    "input": "resources_test/annotation_test_data/TS_Blood_filtered_pseudobulk.h5mu",
    "modality": "rna",
    "output": "output.h5mu",
    "obs_name_filter": "filter_with_counts",
    "obs_count_column": "n_cells",
    "min_count": 15,
    "max_count": 20,
    "output_compression": "gzip",
    "do_subset": False,
}
meta = {"resources_dir": "src/utils/"}
### VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from %s", par["modality"], par["input"])
data = mu.read_h5ad(par["input"], mod=par["modality"])

logger.info("\tUnfiltered data: %s", data)
logger.info("\tComputing aggregations.")


def apply_filter_to_mask(mask, base, filter, comparator):
    new_filt = np.ravel(comparator(base, filter))
    num_removed = np.sum(np.invert(new_filt) & mask)
    mask &= new_filt
    return num_removed, mask


try:
    counts = data.obs[par["obs_count_column"]]
except KeyError:
    raise ValueError(f"Could not find column '{par['obs_count_column']}'")

if not is_numeric_dtype(counts):
    raise ValueError(
        f"Column '{par['obs_count_column']}' does not contain numeric datatype."
    )

if counts.min() < 0:
    raise ValueError(f"Column '{par['obs_count_column']}' contains values < 0.")


# Filter cells
filters = []
if par["min_count"]:
    filters.append(
        (
            "min_count",
            counts,
            ge,
            f"\tRemoving %s observations in {par['obs_count_column']} with <%s counts.",
        )
    )
else:
    logger.info(
        "No minimum count filter applied. Please provide `--min_count` for filtering."
    )

if par["max_count"]:
    filters.append(
        (
            "max_count",
            counts,
            le,
            f"\tRemoving %s observations in {par['obs_count_column']} with >%s counts.",
        )
    )
else:
    logger.info(
        "No maximum count filter applied. Please provide `--max_count` for filtering."
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

if par["do_subset"]:
    modality_data = data[keep_cells, :]
    logger.info("\tFiltered data: %s", modality_data)

logger.info("\tFiltered data: %s", data)
logger.info("Writing output data to %s", par["output"])
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], data, par["output_compression"]
)

logger.info("Finished")
