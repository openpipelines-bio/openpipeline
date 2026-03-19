import sys
import os
import pytest
import subprocess
import re
import mudata as mu
import anndata as ad
import numpy as np
import pandas as pd
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {"resources_dir": "resources_test"}
## VIASH END


@pytest.fixture
def input_mdata(write_mudata_to_file):
    np.random.seed(42)
    n_obs = 100
    cell_types = ["T cell", "B cell", "NK cell", "Monocyte"]
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])

    for method in ["scanvi", "celltypist", "singler"]:
        preds = np.random.choice(cell_types, size=n_obs)
        obs[f"{method}_pred"] = pd.Categorical(preds, categories=cell_types)
        obs[f"{method}_prob"] = np.random.uniform(0.5, 1.0, size=n_obs)

    rna = ad.AnnData(
        X=np.random.rand(n_obs, 10),
        obs=obs,
    )
    prot = ad.AnnData(X=np.random.rand(n_obs, 5), obs=pd.DataFrame(index=obs.index))
    mdata = mu.MuData({"rna": rna, "prot": prot})
    return write_mudata_to_file(mdata)


@pytest.fixture
def tied_mdata(write_mudata_to_file):
    """Fixture where all cells have a guaranteed tie: two methods always disagree."""
    n_obs = 10
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    obs["method_a_pred"] = pd.Categorical(["T cell"] * n_obs)
    obs["method_b_pred"] = pd.Categorical(["B cell"] * n_obs)
    rna = ad.AnnData(X=np.random.rand(n_obs, 5), obs=obs)
    mdata = mu.MuData({"rna": rna})
    return write_mudata_to_file(mdata)


def test_simple_execution(run_component, random_h5mu_path, input_mdata):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_mdata,
            "--input_obs_predictions",
            "scanvi_pred",
            "--input_obs_predictions",
            "celltypist_pred",
            "--input_obs_predictions",
            "singler_pred",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist."

    input_mudata = mu.read_h5mu(input_mdata)
    output_mudata = mu.read_h5mu(output_file)

    assert_annotation_objects_equal(input_mudata.mod["prot"], output_mudata.mod["prot"])

    assert {"consensus_pred", "consensus_score"}.issubset(
        output_mudata.mod["rna"].obs.keys()
    ), "Required keys not found in .obs."

    predictions = output_mudata.mod["rna"].obs["consensus_pred"]
    assert predictions.dtype == "category", (
        "Consensus predictions should be category dtype."
    )
    assert not predictions.isna().all(), "Not all consensus predictions should be NA."

    scores = output_mudata.mod["rna"].obs["consensus_score"]
    assert all(0 <= v <= 1 for v in scores), (
        "Consensus scores should be in the range [0, 1]."
    )


def test_with_weights(run_component, random_h5mu_path, input_mdata):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_mdata,
            "--input_obs_predictions",
            "scanvi_pred",
            "--input_obs_predictions",
            "celltypist_pred",
            "--input_obs_predictions",
            "singler_pred",
            "--weights",
            "1.0",
            "--weights",
            "1.0",
            "--weights",
            "2.0",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist."
    output_mudata = mu.read_h5mu(output_file)

    assert {"consensus_pred", "consensus_score"}.issubset(
        output_mudata.mod["rna"].obs.keys()
    ), "Required keys not found in .obs."


def test_custom_output_obs_names(run_component, random_h5mu_path, input_mdata):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_mdata,
            "--input_obs_predictions",
            "scanvi_pred",
            "--input_obs_predictions",
            "celltypist_pred",
            "--output_obs_predictions",
            "my_consensus_pred",
            "--output_obs_score",
            "my_consensus_score",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist."
    output_mudata = mu.read_h5mu(output_file)

    assert {"my_consensus_pred", "my_consensus_score"}.issubset(
        output_mudata.mod["rna"].obs.keys()
    ), "Custom output keys not found in .obs."


def test_tie_returns_none_by_default(run_component, random_h5mu_path, tied_mdata):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            tied_mdata,
            "--input_obs_predictions",
            "method_a_pred",
            "--input_obs_predictions",
            "method_b_pred",
            "--output",
            output_file,
        ]
    )

    output_mudata = mu.read_h5mu(output_file)
    predictions = output_mudata.mod["rna"].obs["consensus_pred"]
    assert predictions.isna().all(), "All tied predictions should be None."


def test_tie_label(run_component, random_h5mu_path, tied_mdata):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            tied_mdata,
            "--input_obs_predictions",
            "method_a_pred",
            "--input_obs_predictions",
            "method_b_pred",
            "--tie_label",
            "Unknown",
            "--output",
            output_file,
        ]
    )

    output_mudata = mu.read_h5mu(output_file)
    predictions = output_mudata.mod["rna"].obs["consensus_pred"]
    assert (predictions == "Unknown").all(), "All tied predictions should be 'Unknown'."


def test_with_probabilities(run_component, random_h5mu_path, input_mdata):
    output_file = random_h5mu_path()

    run_component(
        [
            "--input",
            input_mdata,
            "--input_obs_predictions",
            "scanvi_pred",
            "--input_obs_predictions",
            "celltypist_pred",
            "--input_obs_predictions",
            "singler_pred",
            "--input_obs_probabilities",
            "scanvi_prob",
            "--input_obs_probabilities",
            "celltypist_prob",
            "--input_obs_probabilities",
            "singler_prob",
            "--output",
            output_file,
        ]
    )

    assert os.path.exists(output_file), "Output file does not exist."
    output_mudata = mu.read_h5mu(output_file)

    assert {"consensus_pred", "consensus_score"}.issubset(
        output_mudata.mod["rna"].obs.keys()
    ), "Required keys not found in .obs."
    assert not output_mudata.mod["rna"].obs["consensus_pred"].isna().all(), (
        "Not all probability-weighted consensus predictions should be NA."
    )


def test_mismatched_probabilities_error(run_component, random_h5mu_path, input_mdata):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_mdata,
                "--input_obs_predictions",
                "scanvi_pred",
                "--input_obs_predictions",
                "celltypist_pred",
                "--input_obs_probabilities",
                "scanvi_prob",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"ValueError: --input_obs_probabilities must have the same length as --input_obs_predictions",
        err.value.stdout.decode("utf-8"),
    )


def test_mismatched_weights_error(run_component, random_h5mu_path, input_mdata):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_mdata,
                "--input_obs_predictions",
                "scanvi_pred",
                "--input_obs_predictions",
                "celltypist_pred",
                "--weights",
                "1.0",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"ValueError: --weights must have the same length as --input_obs_predictions",
        err.value.stdout.decode("utf-8"),
    )


def test_missing_prediction_column_error(run_component, random_h5mu_path, input_mdata):
    output_file = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_mdata,
                "--input_obs_predictions",
                "nonexistent_pred",
                "--output",
                output_file,
            ]
        )
    assert re.search(
        r"ValueError: Prediction column 'nonexistent_pred' not found in .obs",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
