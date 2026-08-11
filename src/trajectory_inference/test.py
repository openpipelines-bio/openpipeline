import sys
import subprocess
import pytest
import mudata as mu
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "name": "trajectory_inference",
    "resources_dir": "resources_test/",
    "executable": "./target/executable/trajectory_inference/trajectory_inference/trajectory_inference",
    "config": "./src/trajectory_inference/config.vsh.yaml",
}
## VIASH END

input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

# This test dataset's leiden clustering was computed on the harmony-integrated
# neighbors graph, and is stored under a non-default column/key name.
OBS_GROUPS = "harmony_integration_leiden_1.0"
NEIGHBORS_KEY = "harmonypy_integration_neighbors"


def test_paga(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--obs_groups",
            OBS_GROUPS,
            "--neighbors_key",
            NEIGHBORS_KEY,
        ]
    )
    assert output_path.is_file()

    input_data = mu.read_h5mu(input_path)
    output_data = mu.read_h5mu(output_path)

    # other modalities should be untouched
    assert_annotation_objects_equal(input_data.mod["prot"], output_data.mod["prot"])

    paga = output_data.mod["rna"].uns["paga"]
    assert "connectivities" in paga
    assert "connectivities_tree" in paga
    assert paga["groups"] == OBS_GROUPS
    n_groups = output_data.mod["rna"].obs[OBS_GROUPS].nunique()
    assert paga["connectivities"].shape == (n_groups, n_groups)


def test_paga_copy(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--obs_groups",
            OBS_GROUPS,
            "--neighbors_key",
            NEIGHBORS_KEY,
            "--copy",
            "True",
        ]
    )
    assert output_path.is_file()

    output_data = mu.read_h5mu(output_path)
    assert "paga" in output_data.mod["rna"].uns


def test_paga_custom_uns_output(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_path,
            "--output",
            output_path,
            "--obs_groups",
            OBS_GROUPS,
            "--neighbors_key",
            NEIGHBORS_KEY,
            "--uns_output",
            "my_paga",
        ]
    )
    assert output_path.is_file()

    output_data = mu.read_h5mu(output_path)
    assert "my_paga" in output_data.mod["rna"].uns
    assert "paga" not in output_data.mod["rna"].uns


def test_paga_obs_groups_not_given_raises(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                "output.h5mu",
                "--neighbors_key",
                NEIGHBORS_KEY,
            ]
        )
    assert (
        "You need to run `tl.leiden` or `tl.louvain` to compute community "
        "labels, or specify `groups='an_existing_key'`"
        in err.value.stdout.decode("utf-8")
    )


def test_paga_invalid_obs_groups_raises(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                "output.h5mu",
                "--obs_groups",
                "does_not_exist",
                "--neighbors_key",
                NEIGHBORS_KEY,
            ]
        )
    assert (
        "ValueError: Requested to use .obs column does_not_exist as the "
        "grouping for PAGA, but the column is not available for modality rna."
        in err.value.stdout.decode("utf-8")
    )


def test_paga_invalid_neighbors_key_raises(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                "output.h5mu",
                "--obs_groups",
                OBS_GROUPS,
                "--neighbors_key",
                "does_not_exist",
            ]
        )
    assert (
        "ValueError: Requested to use .uns key does_not_exist for the "
        "neighbors settings, but the key is not available for modality rna."
        in err.value.stdout.decode("utf-8")
    )


def test_paga_use_rna_velocity_without_graph_raises(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                "output.h5mu",
                "--obs_groups",
                OBS_GROUPS,
                "--neighbors_key",
                NEIGHBORS_KEY,
                "--use_rna_velocity",
                "True",
            ]
        )
    assert (
        "The passed AnnData needs to have an `uns` annotation with key "
        "'velocity_graph' - a sparse matrix from RNA velocity."
        in err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
