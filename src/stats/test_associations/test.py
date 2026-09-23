import sys
import subprocess
import pytest
import numpy as np
import pandas as pd

## VIASH START
meta = {
    "executable": "target/executable/stats/test_associations/test_associations",
    "resources_dir": "src/stats/test_associations/",
    "config": "src/stats/test_associations/config.vsh.yaml",
}
## VIASH END


def _make_tables(tmp_path, n=40, seed=0):
    """One table with a real signal: resp_signal ~ trait_A, resp_null independent."""
    rng = np.random.default_rng(seed)
    trait_a = rng.normal(0, 1, n)
    df = pd.DataFrame(
        {
            "sample_id": [f"s{i:02d}" for i in range(n)],
            "resp_signal": 0.5 * trait_a + rng.normal(0, 0.2, n),
            "resp_null": rng.normal(0, 1, n),
            "trait_A": trait_a,
            "trait_B": rng.normal(0, 1, n),
            "age": rng.integers(60, 90, n).astype(float),
            "sex": rng.choice(["F", "M"], n),
            "cohort": np.where(np.arange(n) < n // 2, "c1", "c2"),
        }
    )
    path = tmp_path / "input.csv"
    df.to_csv(str(path), index=False)
    return df, path


def _make_split_tables(tmp_path, n=40, seed=1):
    """Same data split over a response table and a metadata table."""
    df, _ = _make_tables(tmp_path, n=n, seed=seed)
    responses = df[["sample_id", "resp_signal", "resp_null"]]
    metadata = df[["sample_id", "trait_A", "trait_B", "age", "sex", "cohort"]]
    resp_path = tmp_path / "responses.csv"
    meta_path = tmp_path / "metadata.csv"
    responses.to_csv(str(resp_path), index=False)
    metadata.to_csv(str(meta_path), index=False)
    return resp_path, meta_path


def test_basic(run_component, tmp_path):
    """The true association is found and the null one is not."""
    _, input_path = _make_tables(tmp_path)
    output = tmp_path / "assoc.csv"

    run_component(
        [
            "--input",
            str(input_path),
            "--response_columns",
            "resp_signal",
            "--response_columns",
            "resp_null",
            "--predictor_columns",
            "trait_A",
            "--output",
            str(output),
        ]
    )

    assert output.is_file()
    res = pd.read_csv(str(output))
    assert list(res.columns) == [
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
    assert len(res) == 2, f"expected one row per response, got:\n{res}"
    assert set(res["model"]) == {"ols"}
    assert (res["n"] == 40).all()

    signal = res[res["response"] == "resp_signal"].iloc[0]
    null = res[res["response"] == "resp_null"].iloc[0]
    assert signal["p_value"] < 1e-6, f"signal not detected: {signal.to_dict()}"
    assert null["p_value"] > 0.05, f"null pair came out significant: {null.to_dict()}"
    assert signal["beta"] == pytest.approx(0.5, abs=0.1)
    assert (res["fdr_q"] >= res["p_value"]).all()


def test_join_metadata(run_component, tmp_path):
    """Predictors can live in a second table joined on --join_on."""
    resp_path, meta_path = _make_split_tables(tmp_path)
    output = tmp_path / "assoc_join.csv"

    run_component(
        [
            "--input",
            str(resp_path),
            "--metadata",
            str(meta_path),
            "--join_on",
            "sample_id",
            "--predictor_columns",
            "trait_A",
            "--output",
            str(output),
        ]
    )

    res = pd.read_csv(str(output))
    # --response_columns omitted: the numeric columns of the joined table that are
    # not the predictor, covariate, random-effect or join column
    assert set(res["response"]) == {"resp_signal", "resp_null"}, (
        f"unexpected responses: {sorted(set(res['response']))}"
    )
    assert "trait_B" not in set(res["response"])
    assert (res["n"] == 40).all()


def test_covariates_and_mixed_model(run_component, tmp_path):
    """--random_effect_column switches to mixedlm; covariates are not reported."""
    _, input_path = _make_tables(tmp_path)
    output = tmp_path / "assoc_mixed.csv"

    run_component(
        [
            "--input",
            str(input_path),
            "--response_columns",
            "resp_signal",
            "--predictor_columns",
            "trait_A",
            "--covariate_columns",
            "age",
            "--covariate_columns",
            "sex",
            "--random_effect_column",
            "cohort",
            "--output",
            str(output),
        ]
    )

    res = pd.read_csv(str(output))
    assert set(res["model"]) == {"mixedlm"}
    assert set(res["predictor"]) == {"trait_A"}
    assert not any("age" in t or "sex" in t for t in res["term"]), (
        f"covariates must not be reported: {list(res['term'])}"
    )
    assert res["p_value"].iloc[0] < 1e-6


def test_categorical_predictor(run_component, tmp_path):
    """A string predictor is reported per contrast level."""
    _, input_path = _make_tables(tmp_path)
    output = tmp_path / "assoc_cat.csv"

    run_component(
        [
            "--input",
            str(input_path),
            "--response_columns",
            "resp_signal",
            "--predictor_columns",
            "sex",
            "--output",
            str(output),
        ]
    )

    res = pd.read_csv(str(output))
    assert len(res) == 1
    assert res["term"].iloc[0].startswith("sex["), (
        f"expected a contrast term, got '{res['term'].iloc[0]}'"
    )


def test_transform_logit(run_component, tmp_path):
    """--transform logit changes the fitted coefficients of bounded responses."""
    rng = np.random.default_rng(5)
    n = 40
    trait = rng.normal(0, 1, n)
    prop = 1 / (1 + np.exp(-(0.8 * trait + rng.normal(0, 0.3, n))))
    df = pd.DataFrame({"prop": prop, "trait_A": trait})
    input_path = tmp_path / "prop.csv"
    df.to_csv(str(input_path), index=False)

    betas = {}
    for transform in ("none", "logit"):
        output = tmp_path / f"assoc_{transform}.csv"
        run_component(
            [
                "--input",
                str(input_path),
                "--response_columns",
                "prop",
                "--predictor_columns",
                "trait_A",
                "--transform",
                transform,
                "--output",
                str(output),
            ]
        )
        betas[transform] = pd.read_csv(str(output))["beta"].iloc[0]

    assert betas["none"] != betas["logit"]
    # on the logit scale the generating coefficient (0.8) is recovered
    assert betas["logit"] == pytest.approx(0.8, abs=0.2)


def test_transform_sqrt(run_component, tmp_path):
    """--transform sqrt fits the square root of the response, as BEYOND does."""
    rng = np.random.default_rng(7)
    n = 60
    trait = rng.normal(0, 1, n)
    # Generated linear on the sqrt scale, so only --transform sqrt recovers the slope
    root = 0.4 + 0.1 * trait + rng.normal(0, 0.01, n)
    prop = np.clip(root, 1e-6, None) ** 2
    df = pd.DataFrame({"prop": prop, "trait_A": trait})
    input_path = tmp_path / "prop_sqrt.csv"
    df.to_csv(str(input_path), index=False)

    betas = {}
    for transform in ("none", "sqrt"):
        output = tmp_path / f"assoc_{transform}.csv"
        run_component(
            [
                "--input",
                str(input_path),
                "--response_columns",
                "prop",
                "--predictor_columns",
                "trait_A",
                "--transform",
                transform,
                "--output",
                str(output),
            ]
        )
        betas[transform] = pd.read_csv(str(output))["beta"].iloc[0]

    assert betas["none"] != betas["sqrt"]
    assert betas["sqrt"] == pytest.approx(0.1, abs=0.02), (
        f"sqrt transform did not recover the generating slope: {betas['sqrt']}"
    )


def test_transform_sqrt_negative_response(run_component, tmp_path):
    """--transform sqrt rejects negative responses instead of producing NaN."""
    df = pd.DataFrame(
        {"prop": np.linspace(-0.2, 0.5, 30), "trait_A": np.linspace(0, 1, 30)}
    )
    input_path = tmp_path / "negative.csv"
    df.to_csv(str(input_path), index=False)

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_path),
                "--response_columns",
                "prop",
                "--predictor_columns",
                "trait_A",
                "--transform",
                "sqrt",
                "--output",
                str(tmp_path / "out.csv"),
            ]
        )
    assert "--transform sqrt requires non-negative" in err.value.stdout.decode("utf-8")


def test_fdr_scope(run_component, tmp_path):
    """per_response correction is applied within each response, not over all tests."""
    _, input_path = _make_tables(tmp_path)
    qs = {}
    for scope in ("global", "per_response"):
        output = tmp_path / f"assoc_{scope}.csv"
        run_component(
            [
                "--input",
                str(input_path),
                "--response_columns",
                "resp_signal",
                "--response_columns",
                "resp_null",
                "--predictor_columns",
                "trait_A",
                "--predictor_columns",
                "trait_B",
                "--predictor_columns",
                "age",
                "--fdr_scope",
                scope,
                "--output",
                str(output),
            ]
        )
        res = pd.read_csv(str(output)).set_index(["response", "predictor"])
        qs[scope] = res["fdr_q"]
        assert len(res) == 6

    # 6 tests globally vs 3 per response: the same p-value cannot give the same q
    assert not np.allclose(
        qs["global"].sort_index().to_numpy(), qs["per_response"].sort_index().to_numpy()
    ), "fdr_scope had no effect"


def test_missing_column_raises(run_component, tmp_path):
    """An unknown predictor column is reported with the available columns."""
    _, input_path = _make_tables(tmp_path)
    output = tmp_path / "assoc_err.csv"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_path),
                "--response_columns",
                "resp_signal",
                "--predictor_columns",
                "does_not_exist",
                "--output",
                str(output),
            ]
        )
    assert "not found in the input table" in err.value.stdout.decode("utf-8")


def test_min_observations(run_component, tmp_path):
    """Too few complete observations means no model is fit, and that is an error."""
    df = pd.DataFrame(
        {
            "resp": [0.1, 0.2, 0.3, np.nan, np.nan],
            "trait_A": [1.0, 2.0, 3.0, 4.0, 5.0],
        }
    )
    input_path = tmp_path / "small.csv"
    df.to_csv(str(input_path), index=False)
    output = tmp_path / "assoc_small.csv"

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_path),
                "--response_columns",
                "resp",
                "--predictor_columns",
                "trait_A",
                "--min_observations",
                "5",
                "--output",
                str(output),
            ]
        )
    assert "No models could be fit" in err.value.stdout.decode("utf-8")

    # the same data passes with a lower threshold
    run_component(
        [
            "--input",
            str(input_path),
            "--response_columns",
            "resp",
            "--predictor_columns",
            "trait_A",
            "--min_observations",
            "3",
            "--output",
            str(output),
        ]
    )
    assert pd.read_csv(str(output))["n"].iloc[0] == 3


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
