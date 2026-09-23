import sys
import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
from scipy.stats import f as f_dist

## VIASH START
par = {
    "input": "proportions.csv",
    "pseudotime": "pseudotime.csv",
    "id_column": None,
    "pseudotime_column": "palantir_pseudotime",
    "lam": 0.6,
    "n_pseudotime_bins": 100,
    "min_groups": 5,
    "output": "dynamics.csv",
    "output_stats": None,
}
meta = {
    "resources_dir": "src/trajectory/fit_proportion_dynamics/",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from group_table import read_group_table

logger = setup_logger()


def _read_pseudotime(path, id_column, column):
    """Return a Series of pseudotime indexed by group identifier."""
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"--pseudotime '{path}' has no rows.")
    if id_column not in df.columns:
        raise ValueError(
            f"Identifier column '{id_column}' not found in --pseudotime '{path}'. "
            f"Available: {list(df.columns)}"
        )
    if column not in df.columns:
        raise ValueError(
            f"--pseudotime_column '{column}' not found in --pseudotime '{path}'. "
            f"Available: {list(df.columns)}"
        )
    ids = df[id_column].astype(str)
    if ids.duplicated().any():
        duplicated = sorted(ids[ids.duplicated()].unique())
        raise ValueError(
            f"--pseudotime '{path}' has duplicated identifiers in column "
            f"'{id_column}': {duplicated}"
        )
    values = pd.to_numeric(df[column], errors="coerce")
    series = pd.Series(values.to_numpy(), index=pd.Index(ids.values, name=id_column))
    n_missing = int(series.isna().sum())
    if n_missing:
        logger.warning(
            "Dropping %d group(s) with a missing '%s' value.", n_missing, column
        )
        series = series.dropna()
    return series


def _fit_spline(pt_vals, prop_vals, lam, n_bins):
    """Fit a smoothing spline; return fitted curve, peak, R^2, p-value."""
    x = pt_vals.values if hasattr(pt_vals, "values") else np.asarray(pt_vals)
    y = prop_vals.values if hasattr(prop_vals, "values") else np.asarray(prop_vals)

    # Sort by pseudotime (required by UnivariateSpline)
    order = np.argsort(x)
    x, y = x[order], y[order]

    n = len(x)
    # smoothing_factor s ~ lam * n  (larger = smoother)
    s_factor = lam * n
    spl = UnivariateSpline(x, y, k=3, s=s_factor)

    grid = np.linspace(x.min(), x.max(), n_bins)
    y_fit_grid = spl(grid)
    y_fit_train = spl(x)

    peak_idx = int(np.argmax(y_fit_grid))

    # R^2
    ss_res = float(np.sum((y - y_fit_train) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r_sq = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 0.0

    # F-test: spline vs. intercept-only
    k_spline = max(1, int(spl._data[10]))  # effective degrees of freedom
    df_model = max(1, k_spline - 1)
    df_resid = max(1, n - k_spline)
    ss_mean = ss_tot
    if ss_mean > 1e-12 and df_resid > 0:
        F = ((ss_mean - ss_res) / df_model) / (ss_res / df_resid + 1e-12)
        p_val = float(1.0 - f_dist.cdf(max(0.0, F), df_model, df_resid))
    else:
        p_val = float("nan")

    return {
        "pseudotime_grid": grid,
        "proportion_fitted": y_fit_grid,
        "peak_pseudotime": float(grid[peak_idx]),
        "r_squared": float(np.clip(r_sq, 0.0, 1.0)),
        "p_value": p_val,
    }


def main():
    logger.info("Reading proportions from %s", par["input"])
    prop_df = read_group_table(par["input"], par["id_column"], "--input")
    id_column = prop_df.index.name
    logger.info(
        "Proportions: %d groups x %d labels, identifier column '%s'.",
        prop_df.shape[0],
        prop_df.shape[1],
        id_column,
    )

    logger.info("Reading pseudotime from %s", par["pseudotime"])
    pt = _read_pseudotime(par["pseudotime"], id_column, par["pseudotime_column"])

    common = prop_df.index.intersection(pt.index)
    dropped_prop = prop_df.index.difference(pt.index)
    dropped_pt = pt.index.difference(prop_df.index)
    if len(dropped_prop):
        logger.warning(
            "%d group(s) in --input have no pseudotime and are dropped: %s",
            len(dropped_prop),
            list(dropped_prop[:10]),
        )
    if len(dropped_pt):
        logger.warning(
            "%d group(s) in --pseudotime are absent from --input and are dropped: %s",
            len(dropped_pt),
            list(dropped_pt[:10]),
        )
    if len(common) < par["min_groups"]:
        raise ValueError(
            f"Only {len(common)} group(s) present in both --input and --pseudotime; "
            f"--min_groups is {par['min_groups']}."
        )
    prop_df = prop_df.loc[common]
    pt = pt.loc[common]
    logger.info("Fitting splines for %d labels over %d groups.", prop_df.shape[1], len(common))

    curves = []
    stats = []
    for label in prop_df.columns:
        try:
            result = _fit_spline(
                pt,
                prop_df[label],
                par["lam"],
                par["n_pseudotime_bins"],
            )
        except Exception as exc:
            logger.warning("  Spline fit failed for '%s': %s", label, exc)
            continue
        curves.append(
            pd.DataFrame(
                {
                    "label": label,
                    "pseudotime": result["pseudotime_grid"],
                    "proportion_fitted": result["proportion_fitted"],
                }
            )
        )
        stats.append(
            {
                "label": label,
                "peak_pseudotime": result["peak_pseudotime"],
                "r_squared": result["r_squared"],
                "p_value": result["p_value"],
                "n_groups": len(common),
            }
        )
        logger.info(
            "  %-20s peak_pt=%.3f  R^2=%.3f  p=%.4f",
            label,
            result["peak_pseudotime"],
            result["r_squared"],
            result["p_value"],
        )

    if not curves:
        raise ValueError("No label could be fitted; see the warnings above.")

    curves_df = pd.concat(curves, ignore_index=True)
    curves_df.to_csv(par["output"], index=False)
    logger.info(
        "Written %d fitted curve rows for %d labels to %s.",
        curves_df.shape[0],
        len(stats),
        par["output"],
    )

    if par["output_stats"]:
        pd.DataFrame(stats).to_csv(par["output_stats"], index=False)
        logger.info("Written per-label statistics to %s.", par["output_stats"])


if __name__ == "__main__":
    main()
