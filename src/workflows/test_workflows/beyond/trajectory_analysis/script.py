import sys
import numpy as np
import pytest
from mudata import read_h5mu

##VIASH START
par = {
    "input": "beyond_donors.h5mu",
    "expected_n_obs": 12,
    "expected_n_vars": 9,
    "expected_n_communities": 2,
    "expected_n_terminal_states": 2,
    "expected_n_components": 3,
    "obs_normalize_within": None,
}
meta = {"resources_dir": "resources_test"}
##VIASH END


@pytest.fixture(scope="module")
def adata():
    mdata = read_h5mu(par["input"])
    assert list(mdata.mod) == ["proportions"], (
        f"Expected a single 'proportions' modality, found {list(mdata.mod)}"
    )
    return mdata.mod["proportions"]


def test_one_observation_per_donor(adata):
    # One row per donor - not one per cell, which is the point of the donor-level
    # object: e.g. 152 donors out of 32377 nuclei.
    assert adata.shape == (par["expected_n_obs"], par["expected_n_vars"]), (
        f"Expected {par['expected_n_obs']} donors x {par['expected_n_vars']} "
        f"subpopulations, found {adata.shape}"
    )


def test_proportions(adata):
    X = np.asarray(adata.X)
    np.testing.assert_array_equal(adata.obsm["X_proportions"], X)
    row_sums = X.sum(axis=1)
    within = par["obs_normalize_within"]
    if within is None:
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)
    else:
        # Normalised within the coarser label: a donor's row sums to the number of
        # coarser labels the donor has cells in, not to 1.
        counts = np.asarray(adata.layers["counts"])
        strata = adata.var[within].astype(str).to_numpy()
        n_strata = np.array([len(set(strata[row > 0])) for row in counts], dtype=float)
        np.testing.assert_allclose(row_sums, n_strata, atol=1e-10)
        assert n_strata.max() > 1, "Within-class normalisation had no effect"


def test_landscape_and_pseudotime(adata):
    assert adata.obsm["X_phate"].shape == (
        par["expected_n_obs"],
        par["expected_n_components"],
    )
    pt = adata.obs["palantir_pseudotime"]
    assert pt.notna().all() and pt.between(0, 1).all(), "Pseudotime outside [0, 1]"
    fates = adata.obsm["palantir_fate_probabilities"]
    assert fates.shape[1] == par["expected_n_terminal_states"], (
        f"Expected {par['expected_n_terminal_states']} fate columns, "
        f"found {fates.shape[1]}"
    )


def test_dynamics(adata):
    curves = np.asarray(adata.varm["proportion_dynamics"])
    assert curves.shape[0] == par["expected_n_vars"]
    assert curves.shape[1] == len(adata.uns["proportion_dynamics_pseudotime"])
    for col in ("peak_pseudotime", "r_squared", "p_value", "n_groups"):
        assert f"dynamics_{col}" in adata.var.columns, f"dynamics_{col} missing"


def test_communities(adata):
    communities = adata.var["community_id"]
    assert communities.notna().all()
    assert communities.nunique() == par["expected_n_communities"], (
        f"Expected {par['expected_n_communities']} communities, "
        f"found {communities.nunique()}: {sorted(communities.unique())}"
    )
    for key in ("co_occurrence", "dynamics", "combined"):
        assert adata.varp[f"community_{key}"].shape == (
            par["expected_n_vars"],
            par["expected_n_vars"],
        )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
