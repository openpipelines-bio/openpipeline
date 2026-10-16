import sys
import numpy as np
import pandas as pd
import mudata as mu

## VIASH START
par = {
    "input": "proportions.h5mu",
    "modality": "proportions",
    "response_columns": None,
    "predictor_columns": ["trait_A"],
    "covariate_columns": None,
    "random_effect_column": None,
    "formula": None,
    "transform": "none",
    "pseudocount": 1.0e-6,
    "min_observations": 5,
    "fdr_method": "bh",
    "fdr_scope": "global",
    "output": "associations.csv",
}
meta = {
    "resources_dir": "src/stats/test_associations/",
    "cpus": None,
}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


def _read_input(par):
    """Return (data, responses): one row per observation, responses then .obs."""
    logger.info("Reading %s, modality '%s'.", par["input"], par["modality"])
    adata = mu.read_h5ad(par["input"], mod=par["modality"])
    if adata.n_obs == 0:
        raise ValueError(f"Modality '{par['modality']}' has no observations.")

    var_names = [str(v) for v in adata.var_names]
    responses = list(par["response_columns"] or var_names)
    missing = sorted(set(responses).difference(var_names))
    if missing:
        raise ValueError(
            f"--response_columns {missing} not found in .var_names. "
            f"Available: {var_names[:20]}"
        )

    overlap = sorted(set(responses).intersection(adata.obs.columns))
    if overlap:
        raise ValueError(
            f"Response(s) {overlap} are also .obs column names; rename one of the two."
        )

    idx = [var_names.index(r) for r in responses]
    X = adata.X[:, idx]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    data = pd.DataFrame(X, index=adata.obs_names, columns=responses, dtype=float)
    data = pd.concat([data, adata.obs], axis=1)
    logger.info(
        "%d observations, %d response(s), %d .obs column(s).",
        adata.n_obs,
        len(responses),
        adata.obs.shape[1],
    )
    return data, responses


def _transform_responses(df, responses, method, pseudocount):
    """Apply the requested transformation to the response columns in place."""
    if method == "none":
        return df

    values = df[responses].astype(float)
    if method == "logit":
        if ((values < 0) | (values > 1)).any().any():
            raise ValueError(
                "--transform logit requires all response values in [0, 1]; "
                "found values outside that range."
            )
        clipped = values.clip(lower=pseudocount, upper=1 - pseudocount)
        df[responses] = np.log(clipped / (1 - clipped))
    elif method == "clr":
        if (values < 0).any().any():
            raise ValueError(
                "--transform clr requires non-negative response values; "
                "found negative values."
            )
        shifted = values.where(values > 0, pseudocount)
        log_vals = np.log(shifted)
        df[responses] = log_vals.sub(log_vals.mean(axis=1), axis=0)
    elif method == "sqrt":
        if (values < 0).any().any():
            raise ValueError(
                "--transform sqrt requires non-negative response values; "
                "found negative values."
            )
        df[responses] = np.sqrt(values)
    logger.info(
        "Applied '%s' transform to %d response columns.", method, len(responses)
    )
    return df


def _fit_one(data, response, predictor, par):
    """Fit one model; return a list of result records, one per reported term."""
    import statsmodels.formula.api as smf

    predictor_term = f"Q('{predictor}')"
    if par["formula"]:
        rhs = par["formula"].replace("{predictor}", predictor_term)
    else:
        covariates = par["covariate_columns"] or []
        rhs = " + ".join([predictor_term] + [f"Q('{c}')" for c in covariates])
    formula = f"Q('{response}') ~ {rhs}"

    random_effect = par["random_effect_column"]
    model_name = "mixedlm" if random_effect else "ols"
    try:
        if random_effect:
            fit = smf.mixedlm(formula, data, groups=data[random_effect]).fit(
                reml=True, disp=False
            )
            converged = bool(getattr(fit, "converged", True))
        else:
            fit = smf.ols(formula, data).fit()
            converged = True
    except Exception as exc:
        logger.warning(
            "Model failed for response '%s' x predictor '%s': %s",
            response,
            predictor,
            exc,
        )
        return []

    terms = [t for t in fit.params.index if t.startswith(predictor_term)]
    if not terms:
        logger.warning(
            "No coefficient for predictor '%s' in model of '%s'; skipped.",
            predictor,
            response,
        )
        return []

    records = []
    for term in terms:
        records.append(
            {
                "response": response,
                "predictor": predictor,
                # strip the Q('...') wrapper, keep any categorical level suffix
                "term": term.replace(f"Q('{predictor}')", predictor),
                "beta": float(fit.params[term]),
                "se": float(fit.bse[term]),
                "stat": float(fit.tvalues[term]),
                "p_value": float(fit.pvalues[term]),
                "n": int(len(data)),
                "model": model_name,
                "converged": converged,
            }
        )
    return records


def _correct(results, method, scope):
    """Add an fdr_q column, correcting within the requested scope."""
    if method == "none":
        results["fdr_q"] = results["p_value"]
        return results

    from statsmodels.stats.multitest import multipletests

    method_map = {"bh": "fdr_bh", "bonferroni": "bonferroni"}
    scope_map = {"per_predictor": "predictor", "per_response": "response"}

    results["fdr_q"] = np.nan
    if scope == "global":
        groups = [(None, results.index)]
    else:
        groups = list(results.groupby(scope_map[scope]).groups.items())

    for _, idx in groups:
        _, qvals, _, _ = multipletests(
            results.loc[idx, "p_value"].to_numpy(), method=method_map[method]
        )
        results.loc[idx, "fdr_q"] = qvals

    logger.info("Applied '%s' correction with scope '%s'.", method, scope)
    return results


def main():
    df, responses = _read_input(par)

    predictors = list(par["predictor_columns"])
    covariates = list(par["covariate_columns"] or [])
    random_effect = par["random_effect_column"]

    required = set(predictors) | set(covariates)
    if random_effect:
        required.add(random_effect)
    missing = sorted(c for c in required if c not in df.columns)
    if missing:
        raise ValueError(
            f"Column(s) {missing} not found in .obs. "
            f"Available: {[c for c in df.columns if c not in responses]}"
        )

    df = _transform_responses(df, responses, par["transform"], par["pseudocount"])

    logger.info(
        "Fitting %d response(s) x %d predictor(s).", len(responses), len(predictors)
    )
    records = []
    for response in responses:
        for predictor in predictors:
            model_cols = [response, predictor] + covariates
            if random_effect:
                model_cols.append(random_effect)
            data = df[model_cols].dropna()
            if len(data) < par["min_observations"]:
                logger.debug(
                    "Skipping %s x %s: %d complete observations (min=%d).",
                    response,
                    predictor,
                    len(data),
                    par["min_observations"],
                )
                continue
            records.extend(_fit_one(data, response, predictor, par))

    if not records:
        raise RuntimeError(
            "No models could be fit. Check --min_observations and that the response "
            "and predictor columns share complete observations."
        )

    results = pd.DataFrame(records)
    logger.info("Fit %d model term(s).", len(results))

    results = _correct(results, par["fdr_method"], par["fdr_scope"])
    results = results.sort_values("p_value").reset_index(drop=True)
    results = results[
        [
            "response",
            "predictor",
            "term",
            "beta",
            "se",
            "stat",
            "p_value",
            "fdr_q",
            "n",
            "model",
            "converged",
        ]
    ]

    results.to_csv(par["output"], index=False)
    logger.info("Written %d rows to '%s'.", len(results), par["output"])


if __name__ == "__main__":
    main()
