import sys
import numpy as np
import pandas as pd

## VIASH START
par = {
    "input": "proportions.csv",
    "metadata": "traits.csv",
    "input_join_column": "sample_id",
    "metadata_join_column": None,
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


def _read_table(path, label):
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"{label} '{path}' has no rows.")
    logger.info("%s '%s': %d rows x %d columns.", label, path, *df.shape)
    return df


def _join_tables(df, par):
    """Inner-join the optional metadata table onto the input table."""
    input_col = par["input_join_column"]
    if par["metadata"] is None and input_col is None:
        if par["metadata_join_column"] is not None:
            raise ValueError("--metadata_join_column requires --metadata.")
        return df
    if par["metadata"] is None or input_col is None:
        raise ValueError("--metadata and --input_join_column must be given together.")

    meta_df = _read_table(par["metadata"], "Metadata table")
    meta_col = par["metadata_join_column"] or input_col
    for arg, col, table in (
        ("--input_join_column", input_col, df),
        ("--metadata_join_column", meta_col, meta_df),
    ):
        if col not in table.columns:
            raise ValueError(
                f"{arg} '{col}' not found in its table. "
                f"Available: {list(table.columns)}"
            )
    merged = df.merge(
        meta_df,
        left_on=input_col,
        right_on=meta_col,
        how="inner",
        suffixes=("", "_metadata"),
    )
    if merged.empty:
        raise ValueError(
            f"No rows in common between --input '{input_col}' and "
            f"--metadata '{meta_col}'."
        )
    logger.info(
        "Joined --input '%s' on --metadata '%s': %d rows kept (of %d and %d).",
        input_col,
        meta_col,
        len(merged),
        len(df),
        len(meta_df),
    )
    return merged


def _resolve_responses(df, candidates, par):
    """Return the response columns, defaulting to the unused numeric columns.

    Only columns of --input are candidates: columns that arrive through --metadata
    are predictors, covariates or annotation, never responses.
    """
    if par["response_columns"]:
        return list(par["response_columns"])

    reserved = set(par["predictor_columns"] or [])
    reserved.update(par["covariate_columns"] or [])
    if par["random_effect_column"]:
        reserved.add(par["random_effect_column"])
    if par["input_join_column"]:
        reserved.add(par["input_join_column"])

    responses = [
        c
        for c in candidates
        if c not in reserved and pd.api.types.is_numeric_dtype(df[c])
    ]
    if not responses:
        raise ValueError(
            "--response_columns not given and no unused numeric columns were found "
            f"in --input. Columns: {list(candidates)}"
        )
    logger.info(
        "--response_columns not given; using %d numeric column(s) of --input.",
        len(responses),
    )
    return responses


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
    df = _read_table(par["input"], "Input table")
    input_columns = list(df.columns)
    df = _join_tables(df, par)

    predictors = list(par["predictor_columns"])
    covariates = list(par["covariate_columns"] or [])
    random_effect = par["random_effect_column"]
    responses = _resolve_responses(df, input_columns, par)

    required = set(predictors) | set(covariates) | set(responses)
    if random_effect:
        required.add(random_effect)
    missing = sorted(c for c in required if c not in df.columns)
    if missing:
        raise ValueError(
            f"Column(s) {missing} not found in the input table. "
            f"Available: {list(df.columns)}"
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
