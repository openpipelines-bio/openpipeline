import sys
import numpy as np
import pandas as pd
import mudata as mu
from scipy.interpolate import UnivariateSpline
from scipy.stats import f as f_dist

## VIASH START
par = {
    "input": "test_data/dynamics_input.h5mu",
    "modality": "rna",
    "obs_pseudotime": "palantir_pseudotime",
    "obs_participant_id": "participant_id",
    "uns_proportions": "proportions",
    "n_splines": 8,
    "lam": 0.6,
    "n_pseudotime_bins": 100,
    "min_cells_per_participant": 5,
    "output": "test_data/dynamics_output.h5mu",
    "uns_output": "dynamics",
    "output_compression": None,
}
meta = {
    "resources_dir": "src/trajectory/pseudotime_dynamics/",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def _participant_pseudotimes(adata, participant_col, pseudotime_col, min_cells):
    """Return median pseudotime per participant, dropping those with too few cells."""
    grp = adata.obs.groupby(participant_col, observed=True)[pseudotime_col]
    sizes = grp.count()
    keep = sizes[sizes >= min_cells].index
    excluded = sizes[sizes < min_cells]
    if len(excluded) > 0:
        logger.warning(
            "Excluding %d participant(s) with < %d cells: %s",
            len(excluded),
            min_cells,
            list(excluded.index),
        )
    return grp.median().loc[keep]


def _fit_spline(pt_vals, prop_vals, lam, n_bins):
    """Fit a smoothing spline; return fitted curve, peak, R², p-value."""
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

    # R²
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
        "pseudotime_grid": grid.tolist(),
        "proportion_fitted": y_fit_grid.tolist(),
        "peak_pseudotime": float(grid[peak_idx]),
        "r_squared": float(np.clip(r_sq, 0.0, 1.0)),
        "p_value": p_val,
    }


def main():
    logger.info("Reading input from %s", par["input"])
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]]

    for col in (par["obs_pseudotime"], par["obs_participant_id"]):
        if col not in adata.obs.columns:
            raise ValueError(
                f"Column '{col}' not found in .obs. "
                f"Available: {list(adata.obs.columns)}"
            )
    uns_key = par["uns_proportions"]
    if uns_key not in adata.uns:
        raise ValueError(
            f"Key '{uns_key}' not found in .uns. Available: {list(adata.uns.keys())}"
        )

    # ── participant-level pseudotime ─────────────────────────────────────────
    pt = _participant_pseudotimes(
        adata,
        par["obs_participant_id"],
        par["obs_pseudotime"],
        par["min_cells_per_participant"],
    )
    logger.info("Using %d participants for spline fitting.", len(pt))

    # ── proportion matrix ────────────────────────────────────────────────────
    prop_df = pd.DataFrame(adata.uns[uns_key])  # column-first dict → DataFrame
    common = prop_df.index.intersection(pt.index)
    prop_df = prop_df.loc[common]
    pt = pt.loc[common]

    subpops = prop_df.columns.tolist()
    logger.info("Fitting splines for %d subpopulations.", len(subpops))

    dynamics = {}
    for subpop in subpops:
        try:
            result = _fit_spline(
                pt,
                prop_df[subpop],
                par["lam"],
                par["n_pseudotime_bins"],
            )
            dynamics[subpop] = result
            logger.info(
                "  %-10s peak_pt=%.3f  R²=%.3f  p=%.4f",
                subpop,
                result["peak_pseudotime"],
                result["r_squared"],
                result["p_value"],
            )
        except Exception as exc:
            logger.warning("  Spline fit failed for '%s': %s", subpop, exc)

    adata.uns[par["uns_output"]] = dynamics
    logger.info(
        "Stored dynamics for %d subpopulations in .uns['%s'].",
        len(dynamics),
        par["uns_output"],
    )

    logger.info("Writing output to %s", par["output"])
    write_h5ad_to_h5mu_with_compression(
        output_file=par["output"],
        h5mu=par["input"],
        modality_name=par["modality"],
        modality_data=adata,
        output_compression=par["output_compression"],
    )


if __name__ == "__main__":
    main()
