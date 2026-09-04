import sys
import subprocess
import pytest
import numpy as np
import scipy as sp
import mudata as mu
import h5py
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "name": "paga",
    "resources_dir": "resources_test/",
    "executable": "./target/executable/trajectory_inference/paga/paga",
    "config": "./src/trajectory_inference/paga/config.vsh.yaml",
}
## VIASH END

input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"

# This test dataset's leiden clustering was computed on the harmony-integrated
# neighbors graph, and is stored under a non-default column/key name.
OBS_GROUPS = "harmony_integration_leiden_1.0"
NEIGHBORS_KEY = "harmonypy_integration_neighbors"


def build_paga_args(output, input=input_path, **kwargs):
    """Build a run_component argument list. kwargs map to --key value pairs;
    a kwarg left as None is omitted (used to skip optional CLI flags)."""
    args = ["--input", input, "--output", output]
    for key, value in kwargs.items():
        if value is None:
            continue
        args += [f"--{key}", str(value)]
    return args


def assert_paga_result(output_data, groups_key, input_data=None):
    """Assert that a component run produced a well-formed PAGA result for the
    non-RNA-velocity code path, and return the `paga` dict for further checks."""
    if input_data is not None:
        # other modalities should be untouched
        assert_annotation_objects_equal(input_data.mod["prot"], output_data.mod["prot"])
        # PAGA should only add to .uns -- .obs/.var/.obsm/.X on "rna" itself
        # should otherwise be unaffected (assert_annotation_objects_equal
        # does not compare .uns, so the new "paga" entry doesn't interfere)
        assert_annotation_objects_equal(input_data.mod["rna"], output_data.mod["rna"])

    paga = output_data.mod["rna"].uns["paga"]
    assert "connectivities" in paga
    assert "connectivities_tree" in paga
    assert paga["groups"] == groups_key
    n_groups = output_data.mod["rna"].obs[groups_key].nunique()
    assert paga["connectivities"].shape == (n_groups, n_groups)
    return paga


def assert_component_raises(run_component, args, expected_message):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert expected_message in err.value.stdout.decode("utf-8")


def test_paga(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        build_paga_args(output_path, obs_groups=OBS_GROUPS, uns_neighbors=NEIGHBORS_KEY)
    )
    assert output_path.is_file()

    input_data = mu.read_h5mu(input_path)
    output_data = mu.read_h5mu(output_path)

    assert_paga_result(output_data, OBS_GROUPS, input_data=input_data)


def test_paga_copy(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        build_paga_args(
            output_path,
            obs_groups=OBS_GROUPS,
            uns_neighbors=NEIGHBORS_KEY,
            copy=True,
        )
    )
    assert output_path.is_file()

    output_data = mu.read_h5mu(output_path)
    assert "paga" in output_data.mod["rna"].uns


def test_paga_custom_uns_output(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        build_paga_args(
            output_path,
            obs_groups=OBS_GROUPS,
            uns_neighbors=NEIGHBORS_KEY,
            uns_output="my_paga",
        )
    )
    assert output_path.is_file()

    output_data = mu.read_h5mu(output_path)
    assert "my_paga" in output_data.mod["rna"].uns
    assert "paga" not in output_data.mod["rna"].uns


def test_paga_obs_groups_not_given_raises(run_component):
    assert_component_raises(
        run_component,
        build_paga_args("output.h5mu", uns_neighbors=NEIGHBORS_KEY),
        "--obs_groups is a required argument.",
    )


def test_paga_invalid_obs_groups_raises(run_component):
    assert_component_raises(
        run_component,
        build_paga_args(
            "output.h5mu", obs_groups="does_not_exist", uns_neighbors=NEIGHBORS_KEY
        ),
        "ValueError: Requested to use .obs column does_not_exist as the "
        "grouping for PAGA, but the column is not available for modality rna.",
    )


def test_paga_invalid_uns_neighbors_raises(run_component):
    assert_component_raises(
        run_component,
        build_paga_args(
            "output.h5mu", obs_groups=OBS_GROUPS, uns_neighbors="does_not_exist"
        ),
        "ValueError: Requested to use .uns key does_not_exist for the "
        "neighbors settings, but the key is not available for modality rna.",
    )


def test_paga_use_rna_velocity_without_graph_raises(run_component):
    assert_component_raises(
        run_component,
        build_paga_args(
            "output.h5mu",
            obs_groups=OBS_GROUPS,
            uns_neighbors=NEIGHBORS_KEY,
            use_rna_velocity=True,
        ),
        "The passed AnnData needs to have an `uns` annotation with key "
        "'velocity_graph' - a sparse matrix from RNA velocity.",
    )


def test_paga_v1_0_model(run_component, random_h5mu_path):
    output_file = random_h5mu_path()
    run_component(
        build_paga_args(output_file, obs_groups=OBS_GROUPS, uns_neighbors=NEIGHBORS_KEY)
    )

    output_file_2 = random_h5mu_path()
    run_component(
        build_paga_args(
            output_file_2,
            obs_groups=OBS_GROUPS,
            uns_neighbors=NEIGHBORS_KEY,
            model="v1.0",
        )
    )
    assert output_file.is_file()
    assert output_file_2.is_file()

    paga_mod_v1_2 = assert_paga_result(mu.read_h5mu(output_file), OBS_GROUPS)
    paga_mod_v1_0 = assert_paga_result(mu.read_h5mu(output_file_2), OBS_GROUPS)

    assert not np.array_equal(
        paga_mod_v1_2["connectivities"].toarray(),
        paga_mod_v1_0["connectivities"].toarray(),
    )


def test_paga_rna_velocity(run_component, random_h5mu_path):
    input_data = mu.read_h5mu(input_path)
    n_obs = input_data.mod["rna"].n_obs
    input_data.mod["rna"].uns["velocity_graph"] = sp.sparse.random(
        n_obs, n_obs, density=0.01, format="csr"
    )

    velocity_input_path = random_h5mu_path()
    input_data.write_h5mu(velocity_input_path)

    output_path = random_h5mu_path()
    run_component(
        build_paga_args(
            output_path,
            input=velocity_input_path,
            obs_groups=OBS_GROUPS,
            uns_neighbors=NEIGHBORS_KEY,
            use_rna_velocity=True,
        )
    )
    assert output_path.is_file()

    output_data = mu.read_h5mu(output_path)
    paga = output_data.mod["rna"].uns["paga"]

    # the RNA-velocity branch populates `transitions_confidence` instead of
    # `connectivities` -- this confirms that branch actually ran
    assert "transitions_confidence" in paga
    assert "connectivities" not in paga
    assert paga["groups"] == OBS_GROUPS
    n_groups = output_data.mod["rna"].obs[OBS_GROUPS].nunique()
    assert paga["transitions_confidence"].shape == (n_groups, n_groups)


def test_paga_rna_velocity_custom_uns_key(run_component, random_h5mu_path):
    input_data = mu.read_h5mu(input_path)
    n_obs = input_data.mod["rna"].n_obs
    input_data.mod["rna"].uns["my_velocity_graph"] = sp.sparse.random(
        n_obs, n_obs, density=0.01, format="csr"
    )

    velocity_input_path = random_h5mu_path()
    input_data.write_h5mu(velocity_input_path)

    output_path = random_h5mu_path()
    run_component(
        build_paga_args(
            output_path,
            input=velocity_input_path,
            obs_groups=OBS_GROUPS,
            uns_neighbors=NEIGHBORS_KEY,
            use_rna_velocity=True,
            uns_velocity_graph="my_velocity_graph",
        )
    )
    assert output_path.is_file()

    output_data = mu.read_h5mu(output_path)
    paga = output_data.mod["rna"].uns["paga"]

    # confirms the RNA-velocity branch ran using the graph stored under the
    # custom key, and that the temporary "velocity_graph" alias was cleaned up
    assert "transitions_confidence" in paga
    assert "velocity_graph" not in output_data.mod["rna"].uns


def test_paga_use_rna_velocity_custom_uns_key_missing_raises(run_component):
    assert_component_raises(
        run_component,
        build_paga_args(
            "output.h5mu",
            obs_groups=OBS_GROUPS,
            uns_neighbors=NEIGHBORS_KEY,
            use_rna_velocity=True,
            uns_velocity_graph="does_not_exist",
        ),
        "ValueError: Requested to use .uns key does_not_exist for the RNA "
        "velocity graph, but the key is not available for modality rna.",
    )


def test_paga_default_uns_neighbors(run_component, random_h5mu_path):
    output_default = random_h5mu_path()
    run_component(build_paga_args(output_default, obs_groups=OBS_GROUPS))

    output_key_given = random_h5mu_path()
    run_component(
        build_paga_args(
            output_key_given, obs_groups=OBS_GROUPS, uns_neighbors="neighbors"
        )
    )

    assert output_default.is_file()
    assert output_key_given.is_file()

    paga_default = assert_paga_result(mu.read_h5mu(output_default), OBS_GROUPS)
    paga_key_given = assert_paga_result(mu.read_h5mu(output_key_given), OBS_GROUPS)

    # omitting --uns_neighbors should behave identically to passing the
    # documented default value ("neighbors") explicitly
    assert np.array_equal(
        paga_default["connectivities"].toarray(),
        paga_key_given["connectivities"].toarray(),
    )


def test_paga_modality(run_component, random_h5mu_path):
    input_data = mu.read_h5mu(input_path)
    input_data.mod["rna2"] = input_data.mod["rna"].copy()
    input_data.update()

    input_data_path = random_h5mu_path()
    input_data.write_h5mu(input_data_path)

    output_path = random_h5mu_path()
    run_component(
        build_paga_args(
            output_path,
            input=input_data_path,
            uns_neighbors=NEIGHBORS_KEY,
            obs_groups=OBS_GROUPS,
            modality="rna2",
        )
    )
    assert output_path.is_file()

    output_data = mu.read_h5mu(output_path)

    # the original "rna" modality should be untouched
    assert_annotation_objects_equal(input_data.mod["rna"], output_data.mod["rna"])
    assert "paga" not in output_data.mod["rna"].uns

    paga = output_data.mod["rna2"].uns["paga"]
    assert "connectivities" in paga
    assert "connectivities_tree" in paga
    assert paga["groups"] == OBS_GROUPS
    n_groups = output_data.mod["rna2"].obs[OBS_GROUPS].nunique()
    assert paga["connectivities"].shape == (n_groups, n_groups)


@pytest.mark.parametrize("compression", ["gzip", "lzf"])
def test_paga_output_compression(run_component, random_h5mu_path, compression):
    output_path = random_h5mu_path()
    run_component(
        build_paga_args(
            output_path,
            obs_groups=OBS_GROUPS,
            uns_neighbors=NEIGHBORS_KEY,
            output_compression=compression,
        )
    )
    assert output_path.is_file()

    # content must still be correct after the compress-and-rewrite step
    output_data = mu.read_h5mu(output_path)
    assert_paga_result(output_data, OBS_GROUPS)

    # every non-scalar dataset newly written by this component must carry the
    # requested codec. Pre-existing datasets from the input file (e.g. under
    # mod/prot) may already have their own compression, which compress_h5mu
    # intentionally leaves untouched, so those are out of scope here.
    with h5py.File(output_path, "r") as f:

        def check_compression(name, obj):
            if not name.startswith("mod/rna/uns/paga"):
                return
            if isinstance(obj, h5py.Dataset) and obj.shape not in (None, ()):
                assert obj.compression == compression, (
                    f"{name} is not {compression}-compressed"
                )

        f.visititems(check_compression)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
