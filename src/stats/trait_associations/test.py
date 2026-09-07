import sys
import pytest
import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu

## VIASH START
meta = {
    "executable": "target/executable/stats/trait_associations/trait_associations",
    "resources_dir": "src/stats/trait_associations/",
    "config": "src/stats/trait_associations/config.vsh.yaml",
}
## VIASH END


def _make_input(tmp_path, seed=0, n_participants=30, n_subpops=6):
    """Create a synthetic h5mu and traits CSV for testing."""
    rng = np.random.default_rng(seed)
    subpops = [f"SP{i}" for i in range(n_subpops)]
    participant_ids = [f"donor_{i:02d}" for i in range(n_participants)]

    # Proportions: each row sums to 1
    raw = rng.dirichlet(np.ones(n_subpops), size=n_participants)
    prop_df = pd.DataFrame(raw, index=participant_ids, columns=subpops)

    # Traits: trait_A correlated with SP0, trait_B random, cohort for random effect
    trait_a = prop_df["SP0"].values + rng.normal(0, 0.05, n_participants)
    trait_b = rng.normal(0, 1, n_participants)
    cohort = np.array(
        ["cohort_A"] * (n_participants // 2)
        + ["cohort_B"] * (n_participants - n_participants // 2)
    )

    traits_csv = tmp_path / "traits.csv"
    traits_df = pd.DataFrame(
        {
            "participant_id": participant_ids,
            "trait_A": trait_a,
            "trait_B": trait_b,
            "age": rng.integers(60, 90, n_participants).astype(float),
            "cohort": cohort,
        }
    )
    traits_df.to_csv(str(traits_csv), index=False)

    # Build AnnData
    n_cells = n_participants * 5
    obs_rows = []
    for i, pid in enumerate(participant_ids):
        sp = rng.choice(subpops, 5, p=raw[i]).tolist()
        for s in sp:
            obs_rows.append({"participant_id": pid, "subpopulation": s})
    obs = pd.DataFrame(obs_rows, index=[f"c{i}" for i in range(n_cells)])
    adata = ad.AnnData(
        X=rng.integers(0, 50, (n_cells, 3)).astype("float32"),
        obs=obs,
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
    )
    adata.uns["proportions"] = prop_df.to_dict()

    h5mu_path = tmp_path / "input.h5mu"
    mu.MuData({"rna": adata}).write_h5mu(str(h5mu_path))
    return h5mu_path, traits_csv


def test_basic_association(run_component, tmp_path):
    """Basic run: two traits, outputs stored in uns and associations detected."""
    h5mu, traits_csv = _make_input(tmp_path)
    out = tmp_path / "out.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--traits_csv",
            str(traits_csv),
            "--trait_columns",
            "trait_A",
            "--trait_columns",
            "trait_B",
            "--output",
            str(out),
        ]
    )

    assert out.is_file()
    adata = mu.read_h5mu(str(out)).mod["rna"]
    assert "trait_associations" in adata.uns
    res = pd.DataFrame(adata.uns["trait_associations"])
    assert "subpopulation" in res.columns
    assert "trait" in res.columns
    assert "beta" in res.columns
    assert "p_value" in res.columns
    assert "fdr_q" in res.columns
    assert len(res) > 0


def test_trait_a_detected(run_component, tmp_path):
    """trait_A (correlated with SP0) should show significant association for SP0."""
    h5mu, traits_csv = _make_input(tmp_path, seed=42)
    out = tmp_path / "out2.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--traits_csv",
            str(traits_csv),
            "--trait_columns",
            "trait_A",
            "--output",
            str(out),
        ]
    )

    res = pd.DataFrame(mu.read_h5mu(str(out)).mod["rna"].uns["trait_associations"])
    sp0_row = res[res["subpopulation"] == "SP0"]
    assert len(sp0_row) == 1
    # SP0 was generated to correlate with trait_A -> beta should be positive and p small
    assert float(sp0_row["beta"].iloc[0]) > 0
    assert float(sp0_row["p_value"].iloc[0]) < 0.05


def test_random_effect(run_component, tmp_path):
    """Mixed model with random effect runs without error."""
    h5mu, traits_csv = _make_input(tmp_path, seed=1)
    out = tmp_path / "out3.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--traits_csv",
            str(traits_csv),
            "--trait_columns",
            "trait_A",
            "--random_effect_column",
            "cohort",
            "--output",
            str(out),
        ]
    )

    adata = mu.read_h5mu(str(out)).mod["rna"]
    assert "trait_associations" in adata.uns
    res = pd.DataFrame(adata.uns["trait_associations"])
    assert len(res) > 0


def test_covariates(run_component, tmp_path):
    """Covariates are accepted and model still runs."""
    h5mu, traits_csv = _make_input(tmp_path, seed=2)
    out = tmp_path / "out4.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--traits_csv",
            str(traits_csv),
            "--trait_columns",
            "trait_A",
            "--covariate_columns",
            "age",
            "--output",
            str(out),
        ]
    )

    res = pd.DataFrame(mu.read_h5mu(str(out)).mod["rna"].uns["trait_associations"])
    assert len(res) > 0


def test_output_csv(run_component, tmp_path):
    """CSV output is written when --output_csv is specified."""
    h5mu, traits_csv = _make_input(tmp_path, seed=3)
    out = tmp_path / "out5.h5mu"
    csv_out = tmp_path / "results.csv"

    run_component(
        [
            "--input",
            str(h5mu),
            "--traits_csv",
            str(traits_csv),
            "--trait_columns",
            "trait_A",
            "--output",
            str(out),
            "--output_csv",
            str(csv_out),
        ]
    )

    assert csv_out.is_file()
    df = pd.read_csv(str(csv_out))
    assert "subpopulation" in df.columns
    assert "p_value" in df.columns


def test_fdr_bonferroni(run_component, tmp_path):
    """Bonferroni correction runs without error."""
    h5mu, traits_csv = _make_input(tmp_path, seed=4)
    out = tmp_path / "out6.h5mu"

    run_component(
        [
            "--input",
            str(h5mu),
            "--traits_csv",
            str(traits_csv),
            "--trait_columns",
            "trait_A",
            "--trait_columns",
            "trait_B",
            "--fdr_method",
            "bonferroni",
            "--output",
            str(out),
        ]
    )

    res = pd.DataFrame(mu.read_h5mu(str(out)).mod["rna"].uns["trait_associations"])
    assert len(res) > 0


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
