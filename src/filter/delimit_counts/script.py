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
    "obs_count_column": ["n_cells"],
    "obs_name_filter": ["filter_with_counts"],
    "var_count_column": None,
    "var_name_filter": None,
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


def delimit_axis(dataframe, axis_name, count_columns, name_filters):
    """Threshold one or more columns of a dataframe (.obs or .var) and store the
    resulting boolean masks. Returns the combined (AND) mask across all columns."""
    count_columns = count_columns or []
    name_filters = name_filters or []
    if len(count_columns) != len(name_filters):
        raise ValueError(
            f"The number of --{axis_name}_count_column values ({len(count_columns)}) "
            f"must match the number of --{axis_name}_name_filter values ({len(name_filters)})."
        )

    combined_mask = np.repeat(True, dataframe.shape[0])
    for count_column, name_filter in zip(count_columns, name_filters):
        try:
            counts = dataframe[count_column]
        except KeyError:
            raise ValueError(f"Could not find column '{count_column}' in .{axis_name}")

        if not is_numeric_dtype(counts):
            raise ValueError(
                f"Column '{count_column}' does not contain numeric datatype."
            )
        if counts.min() < 0:
            raise ValueError(f"Column '{count_column}' contains values < 0.")

        filters = []
        if par["min_count"] is not None:
            filters.append(
                (
                    par["min_count"],
                    ge,
                    f"\tRemoving %s entries in {count_column} with <%s counts.",
                )
            )
        if par["max_count"] is not None:
            filters.append(
                (
                    par["max_count"],
                    le,
                    f"\tRemoving %s entries in {count_column} with >%s counts.",
                )
            )
        if not filters:
            logger.info(
                "No filter applied to column '%s'. Please provide `--min_count` "
                "and/or `--max_count` for filtering.",
                count_column,
            )

        keep = np.repeat(True, dataframe.shape[0])
        for threshold, comparator, message in filters:
            num_removed, keep = apply_filter_to_mask(
                keep, counts, threshold, comparator
            )
            logger.info(message, num_removed, threshold)

        dataframe[name_filter] = keep
        combined_mask &= keep

    return combined_mask


if not par["obs_count_column"] and not par["var_count_column"]:
    raise ValueError(
        "At least one of --obs_count_column or --var_count_column must be provided."
    )

keep_obs = delimit_axis(
    data.obs, "obs", par["obs_count_column"], par["obs_name_filter"]
)
keep_var = delimit_axis(
    data.var, "var", par["var_count_column"], par["var_name_filter"]
)

if par["do_subset"]:
    data = data[keep_obs, keep_var]
    logger.info("\tFiltered data: %s", data)

logger.info("Writing output data to %s", par["output"])
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], data, par["output_compression"]
)

logger.info("Finished")
