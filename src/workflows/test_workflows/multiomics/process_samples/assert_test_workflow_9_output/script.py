import sys
import mudata as mu
import pytest

##VIASH START
par = {
    "input": "output.h5mu",
    "orig_input": [
        "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    ],
}

meta = {"resources_dir": "src/base"}
##VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger

logger = setup_logger()


@pytest.fixture
def input_h5mu():
    return mu.read_h5mu(par["orig_input"][0])


@pytest.fixture
def output_h5mu():
    return mu.read_h5mu(par["input"])


def test_modalities_present(output_h5mu):
    assert all(k in output_h5mu.mod.keys() for k in ("rna", "prot"))


def test_intersect_obs_applied(output_h5mu):
    rna_obs = set(output_h5mu["rna"].obs_names)
    prot_obs = set(output_h5mu["prot"].obs_names)
    assert rna_obs == prot_obs, (
        "After intersect_obs, RNA and prot modalities should share the same observations. "
        f"RNA-only: {len(rna_obs - prot_obs)}, prot-only: {len(prot_obs - rna_obs)}."
    )


def test_intersection_is_subset_of_input(output_h5mu, input_h5mu):
    input_rna_obs = set(input_h5mu["rna"].obs_names)
    input_prot_obs = set(input_h5mu["prot"].obs_names)
    expected_intersection = input_rna_obs & input_prot_obs
    # observations were renamed by add_id_make_observation_keys_unique ("pbmc_<barcode>")
    output_rna_obs = {
        name.removeprefix("pbmc_") for name in output_h5mu["rna"].obs_names
    }
    assert output_rna_obs.issubset(expected_intersection), (
        "Output observations should be a subset of the intersection of input modality observations."
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "--import-mode=importlib"]))
