import sys
import mudata as mu
import pytest
import os

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

# test_wf10 applies per-modality cell filtering (rna_min_counts=1000,
# prot_min_counts=100) without enabling --intersect_obs. The checks below
# also guard against intersect_obs being silently applied in process_batches:
# if it were, rna and prot would end up with the same observation set.
MODALITIES = ("rna", "prot")


@pytest.fixture
def input_h5mu():
    return mu.read_h5mu(par["orig_input"][0])


@pytest.fixture
def output_h5mu():
    return mu.read_h5mu(par["input"])


def test_filtering_removes_cells(output_h5mu, input_h5mu):
    for modality in MODALITIES:
        assert output_h5mu[modality].n_obs < input_h5mu[modality].n_obs, (
            f"Modality {modality!r}: expected filtering to reduce n_obs "
            f"below the input ({input_h5mu[modality].n_obs}), but got "
            f"{output_h5mu[modality].n_obs}."
        )


def test_modalities_have_different_obs_counts(output_h5mu):
    # Regression guard: with independent per-modality thresholds and
    # --intersect_obs disabled, rna and prot must not share n_obs. Equal
    # counts would indicate intersect_obs was silently applied.
    rna_n = output_h5mu["rna"].n_obs
    prot_n = output_h5mu["prot"].n_obs
    assert rna_n != prot_n, (
        f"rna and prot modalities have the same n_obs ({rna_n}); "
        f"this suggests intersect_obs was silently applied."
    )


if __name__ == "__main__":
    sys.exit(
        pytest.main(
            [
                "--import-mode=importlib",
                "-o",
                "python_files=script*.py .viash_script.py",
                os.path.dirname(__file__),
            ]
        )
    )
