import sys
import numpy as np
import pandas as pd
import mudata as mu

## VIASH START
par = {
    "input": "input.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obs_column": "random_normal_int",
    "obs_min_quantile": 0.1,
    "obs_max_quantile": 0.9,
    "var_column": "gene_symbol",
    "var_min_quantile": 0.05,
    "var_max_quantile": 0.95,
    "obs_name_filter": "filter_with_quantile",
    "var_name_filter": "filter_with_quantile",
    "do_subset": False,
    "min_values_for_quantile": 10,
    "obs_log1p_transform": False,
    "var_log1p_transform": False,
    "obs_log1p_column": None,
    "var_log1p_column": None,
}

meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def filter_by_quantile(
    values, min_quantile=None, max_quantile=None, min_values_for_quantile=10
):
    """
    Create a boolean mask based on quantile thresholds.

    Parameters:
    -----------
    values : array-like
        Numeric values to filter
    min_quantile : float, optional
        Minimum quantile threshold (0.0-1.0). Values below this quantile are filtered out.
    max_quantile : float, optional
        Maximum quantile threshold (0.0-1.0). Values above this quantile are filtered out.
    min_values_for_quantile : int, default 10
        Minimum number of values required to apply quantile filtering. If fewer values are present,
        filtering will be skipped to prevent data loss.

    Returns:
    --------
    mask : array
        Boolean array where True indicates the value should be kept
    """
    mask = np.ones(len(values), dtype=bool)
    n_values = len(values)

    # Check for NaN values and raise error if found
    if pd.isna(values).any():
        raise ValueError(
            "Column contains NaN values. Please clean the data before applying quantile filtering."
        )

    # Check if there's sufficient data for quantile filtering
    if n_values < min_values_for_quantile:
        logger.info(
            f"Insufficient data ({n_values} values) for quantile filtering. "
            f"Minimum required: {min_values_for_quantile}. Keeping all data."
        )
        return mask

    if min_quantile is not None:
        logger.info(f"Applying minimum quantile filter: {min_quantile}")
        threshold = np.quantile(values, min_quantile)
        mask &= values >= threshold
        logger.info(
            f"{(values < threshold).sum()}/{len(mask)} values are below the minimum quantile threshold"
        )

    if max_quantile is not None:
        logger.info(f"Applying maximum quantile filter: {max_quantile}")
        threshold = np.quantile(values, max_quantile)
        mask &= values <= threshold
        logger.info(
            f"{(values > threshold).sum()}/{len(mask)} values are above the maximum quantile threshold"
        )

    filtered_count = np.sum(~mask)
    logger.info(f"{filtered_count}/{len(mask)} values will be removed")

    return mask


def main():
    logger = setup_logger()

    logger.info(f"Loading input file: {par['input']}")
    adata = mu.read_h5ad(par["input"], mod=par["modality"])

    filter_config = {
        "obs": {
            "filter": "observations",
            "column": par["obs_column"],
            "min_quantile": par["obs_min_quantile"],
            "max_quantile": par["obs_max_quantile"],
            "name_filter": par["obs_name_filter"],
            "log1p_transform": par["obs_log1p_transform"],
            "log1p_column": par["obs_log1p_column"],
        },
        "var": {
            "filter": "features",
            "column": par["var_column"],
            "min_quantile": par["var_min_quantile"],
            "max_quantile": par["var_max_quantile"],
            "name_filter": par["var_name_filter"],
            "log1p_transform": par["var_log1p_transform"],
            "log1p_column": par["var_log1p_column"],
        },
    }

    filter_results = {}
    for filter_type, config in filter_config.items():
        if config["column"] is not None:
            logger.info(
                f"Applying {config['filter']} filtering on column: {config['column']}"
            )

            data_frame = getattr(adata, filter_type)
            if config["column"] not in data_frame.columns:
                raise ValueError(
                    f"Column '{config['column']}' not found in .{filter_type}. Available columns: {list(data_frame.columns)}"
                )

            if config["min_quantile"] is None and config["max_quantile"] is None:
                logger.warning(
                    f"No `--{filter_type}_min_quantile` or `--{filter_type}_max_quantile` provided, no {config['filter']} filtering will be applied."
                )
                filtered_values = np.ones(data_frame.shape[0], dtype=bool)  # Keep all
            else:
                values = data_frame[config["column"]].values

                # Check if column contains numeric data
                if not pd.api.types.is_numeric_dtype(values):
                    raise ValueError(
                        f"Column '{config['column']}' must contain numeric data for quantile filtering"
                    )

                # Apply log1p transformation if requested
                if config["log1p_transform"]:
                    logger.info(f"Applying log1p transformation to {config['column']}")
                    # Use custom column name if provided, otherwise default to log1p_{column}
                    log1p_column_name = (
                        config.get("log1p_column") or f"log1p_{config['column']}"
                    )
                    log1p_values = np.log1p(values)
                    data_frame[log1p_column_name] = log1p_values
                    logger.info(
                        f"Stored log1p transformed values in .{filter_type}['{log1p_column_name}']"
                    )
                    # Use transformed values for quantile filtering
                    filter_values = log1p_values
                else:
                    # Use original values for quantile filtering
                    filter_values = values

                filtered_values = filter_by_quantile(
                    filter_values,
                    min_quantile=config["min_quantile"],
                    max_quantile=config["max_quantile"],
                    min_values_for_quantile=par["min_values_for_quantile"],
                )

            # Store filter masks
            logger.info(
                f"Storing {config['filter']} filter in .{filter_type}['{config['name_filter']}']"
            )
            data_frame[config["name_filter"]] = filtered_values
            filter_results[f"{filter_type}_filter"] = filtered_values

    # Apply subsetting if requested
    if par["do_subset"]:
        logger.info("Subsetting data based on filter masks")
        original_shape = adata.shape
        obs_mask = filter_results.get("obs_filter", np.ones(adata.n_obs, dtype=bool))
        var_mask = filter_results.get("var_filter", np.ones(adata.n_vars, dtype=bool))
        adata = adata[obs_mask, var_mask].copy()
        logger.info(f"Data shape: {original_shape} -> {adata.shape}")

    # Write output
    logger.info(f"Writing output to: {par['output']}")
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], adata, par["output_compression"]
    )
    logger.info("Filtering complete!")


if __name__ == "__main__":
    main()
