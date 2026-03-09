#!/usr/bin/env python3

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
}

meta = {"resources_dir": "src/utils"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def filter_by_quantile(values, min_quantile=None, max_quantile=None):
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

    Returns:
    --------
    mask : array
        Boolean array where True indicates the value should be kept
    """
    mask = np.ones(len(values), dtype=bool)

    # Check for NaN values and raise error if found
    if pd.isna(values).any():
        raise ValueError(
            "Column contains NaN values. Please clean the data before applying quantile filtering."
        )

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
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]].copy()

    # Apply obs filtering
    if par["obs_column"] is not None:
        logger.info(f"Applying observation filtering on column: {par['obs_column']}")

        if par["obs_column"] not in adata.obs.columns:
            raise ValueError(
                f"Column '{par['obs_column']}' not found in .obs. Available columns: {list(adata.obs.columns)}"
            )

        if par["obs_min_quantile"] is None and par["obs_max_quantile"] is None:
            logger.warning(
                "No `--obs_min_quantile` or `--obs_max_quantile` provided, so .obs filtering will be applied."
            )
            obs_filter = np.ones(adata.n_obs, dtype=bool)  # Keep all observations

        else:
            values = adata.obs[par["obs_column"]].values

            # Check if column contains numeric data
            if not pd.api.types.is_numeric_dtype(values):
                raise ValueError(
                    f"Column '{par['obs_column']}' must contain numeric data for quantile filtering"
                )

            obs_filter = filter_by_quantile(
                values,
                min_quantile=par["obs_min_quantile"],
                max_quantile=par["obs_max_quantile"],
            )

        # Store filter masks
        logger.info(f"Storing observation filter in .obs['{par['obs_name_filter']}']")
        adata.obs[par["obs_name_filter"]] = obs_filter

    # Apply var filtering
    if par["var_column"] is not None:
        logger.info(f"Applying variable filtering on column: {par['var_column']}")

        if par["var_column"] not in adata.var.columns:
            raise ValueError(
                f"Column '{par['var_column']}' not found in .var. Available columns: {list(adata.var.columns)}"
            )

        if par["var_min_quantile"] is None and par["var_max_quantile"] is None:
            logger.warning(
                "No `--var_min_quantile` or `--var_max_quantile` provided, no .var filtering will be applied."
            )
            var_filter = np.ones(adata.n_vars, dtype=bool)  # Keep all features

        else:
            values = adata.var[par["var_column"]].values

            # Check if column contains numeric data
            if not pd.api.types.is_numeric_dtype(values):
                raise ValueError(
                    f"Column '{par['var_column']}' must contain numeric data for quantile filtering"
                )

            var_filter = filter_by_quantile(
                values,
                min_quantile=par["var_min_quantile"],
                max_quantile=par["var_max_quantile"],
            )

        # Store filter masks
        logger.info(f"Storing variable filter in .var['{par['var_name_filter']}']")
        adata.var[par["var_name_filter"]] = var_filter

    # Apply subsetting if requested
    if par["do_subset"]:
        logger.info("Subsetting data based on filter masks")
        original_shape = adata.shape
        adata = adata[obs_filter, var_filter].copy()
        logger.info(f"Data shape: {original_shape} -> {adata.shape}")

    # Update the modality in mdata
    mdata.mod[par["modality"]] = adata

    # Write output
    logger.info(f"Writing output to: {par['output']}")
    mdata.write_h5mu(par["output"])

    logger.info("Filtering complete!")


if __name__ == "__main__":
    main()
