import sys
import numpy as np
import pandas as pd
import mudata as mu

## VIASH START
par = {
    "input": "test_data/traits_input.h5mu",
    "modality": "rna",
    "uns_proportions": "proportions",
    "traits_csv": "test_data/traits.csv",
    "participant_id_column": "participant_id",
    "trait_columns": ["trait_A", "trait_B"],
    "covariate_columns": None,
    "random_effect_column": None,
    "fdr_method": "bh",
    "min_participants": 5,
    "output": "test_data/traits_output.h5mu",
    "output_csv": None,
    "uns_output": "trait_associations",
    "output_compression": None,
}
meta = {
    "resources_dir": "src/stats/trait_associations/",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()


def _fit_model(sub_df, trait_col, covariate_cols, random_effect_col):
    """Fit one (subpopulation, trait) association model.

    Returns dict with beta, se, t_stat, p_value, n_participants, converged.
    Returns None when the model cannot be fit.
    """
    import statsmodels.formula.api as smf

    data = sub_df[["proportion", trait_col] + covariate_cols].dropna()
    if len(data) < 2:
        return None

    # Build RHS of formula (trait + covariates)
    rhs_terms = [f"Q('{trait_col}')"] + [f"Q('{c}')" for c in covariate_cols]
    formula = "proportion ~ " + " + ".join(rhs_terms)

    try:
        if random_effect_col and random_effect_col in sub_df.columns:
            data = sub_df[["proportion", trait_col, random_effect_col] + covariate_cols].dropna()
            if len(data) < 2:
                return None
            result = smf.mixedlm(formula, data, groups=data[random_effect_col]).fit(
                reml=True, disp=False
            )
        else:
            result = smf.ols(formula, data).fit()

        trait_key = f"Q('{trait_col}')"
        params = result.params
        bse = result.bse
        tvalues = result.tvalues
        pvalues = result.pvalues

        if trait_key not in params.index:
            return None

        return {
            "beta": float(params[trait_key]),
            "se": float(bse[trait_key]),
            "t_stat": float(tvalues[trait_key]),
            "p_value": float(pvalues[trait_key]),
            "n_participants": int(len(data)),
            "converged": True,
        }
    except Exception as exc:
        logger.debug("Model failed for trait '%s': %s", trait_col, exc)
        return None


def main():
    logger.info("Reading input from %s", par["input"])
    mdata = mu.read_h5mu(par["input"])
    adata = mdata.mod[par["modality"]]

    # ── validate ──────────────────────────────────────────────────────────────
    uns_key = par["uns_proportions"]
    if uns_key not in adata.uns:
        raise ValueError(
            f"Key '{uns_key}' not found in .uns. Available: {list(adata.uns.keys())}"
        )

    # ── load proportion matrix ────────────────────────────────────────────────
    prop_df = pd.DataFrame(adata.uns[uns_key])  # participants × subpopulations
    logger.info(
        "Proportion matrix: %d participants × %d subpopulations.",
        *prop_df.shape,
    )

    # ── load traits table ─────────────────────────────────────────────────────
    traits_df = pd.read_csv(par["traits_csv"])
    pid_col = par["participant_id_column"]
    if pid_col not in traits_df.columns:
        raise ValueError(
            f"Participant ID column '{pid_col}' not found in traits CSV. "
            f"Available: {list(traits_df.columns)}"
        )
    traits_df = traits_df.set_index(pid_col)

    trait_cols = par["trait_columns"] or []
    for tc in trait_cols:
        if tc not in traits_df.columns:
            raise ValueError(
                f"Trait column '{tc}' not found in traits CSV. "
                f"Available: {list(traits_df.columns)}"
            )

    covariate_cols = par["covariate_columns"] or []
    random_effect_col = par["random_effect_column"]

    # ── merge proportions with traits ─────────────────────────────────────────
    common_participants = prop_df.index.intersection(traits_df.index)
    if len(common_participants) == 0:
        raise ValueError(
            "No participants found in common between proportion matrix index "
            f"and '{pid_col}' column in traits CSV."
        )
    logger.info(
        "%d participants in common between proportion matrix and traits CSV.",
        len(common_participants),
    )

    prop_sub = prop_df.loc[common_participants]
    traits_sub = traits_df.loc[common_participants]

    # ── run association models ────────────────────────────────────────────────
    records = []
    subpops = prop_sub.columns.tolist()
    min_n = par["min_participants"]

    for subpop in subpops:
        for trait in trait_cols:
            merged = pd.concat(
                [prop_sub[[subpop]].rename(columns={subpop: "proportion"}), traits_sub],
                axis=1,
            )
            n_complete = merged[[trait]].dropna().__len__()
            if n_complete < min_n:
                logger.debug(
                    "Skipping %s × %s: only %d complete cases (min=%d).",
                    subpop, trait, n_complete, min_n,
                )
                continue

            fit = _fit_model(merged, trait, covariate_cols, random_effect_col)
            if fit is None:
                logger.warning("Model failed for subpopulation '%s' × trait '%s'.", subpop, trait)
                continue

            records.append({
                "subpopulation": subpop,
                "trait": trait,
                **fit,
            })

    if not records:
        raise RuntimeError(
            "No association models could be fit. Check that trait and proportion "
            "matrix share participants, and that --min_participants is not too high."
        )

    results_df = pd.DataFrame(records)
    logger.info("Fit %d association models.", len(results_df))

    # ── FDR correction ────────────────────────────────────────────────────────
    fdr_method = par["fdr_method"]
    if fdr_method == "none":
        results_df["fdr_q"] = results_df["p_value"]
    else:
        from statsmodels.stats.multitest import multipletests
        method_map = {"bh": "fdr_bh", "bonferroni": "bonferroni"}
        _, fdr_q, _, _ = multipletests(
            results_df["p_value"].values, method=method_map[fdr_method]
        )
        results_df["fdr_q"] = fdr_q
        logger.info("Applied '%s' FDR correction.", fdr_method)

    results_df = results_df.sort_values("p_value").reset_index(drop=True)

    # ── store in uns ──────────────────────────────────────────────────────────
    adata.uns[par["uns_output"]] = results_df.to_dict(orient="list")
    logger.info(
        "Stored %d association results in .uns['%s'].",
        len(results_df),
        par["uns_output"],
    )

    # ── optional CSV output ───────────────────────────────────────────────────
    if par.get("output_csv"):
        results_df.to_csv(par["output_csv"], index=False)
        logger.info("Written CSV to %s", par["output_csv"])

    # ── write output ──────────────────────────────────────────────────────────
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
