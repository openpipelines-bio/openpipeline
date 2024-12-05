import pandas as pd
import numpy as np
import mudata as mu
import anndata as ad
import re
import sys
import pytest
import subprocess


@pytest.fixture
def input_h5mu():
    # generate data
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]], index=df.columns, columns=["Feat"])
    ad1 = ad.AnnData(df, obs=obs, var=var)
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = ad.AnnData(df, obs=obs2, var=var2)
    tmp_mudata = mu.MuData({"mod1": ad1, "mod2": ad2})
    return tmp_mudata


@pytest.fixture
def input_h5mu_path(write_mudata_to_file, input_h5mu):
    return write_mudata_to_file(input_h5mu)


def test_copy_obs(run_component, random_h5mu_path, input_h5mu, input_h5mu_path):
    output_h5mu_path = random_h5mu_path()

    args = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--modality",
        "mod1",
        "--input_obs_key",
        "Obs",
        "--output_obs_key",
        "Obs_copy",
    ]

    run_component(args)

    assert output_h5mu_path.is_file()

    output_h5mu = mu.read_h5mu(output_h5mu_path)

    assert "Obs_copy" in output_h5mu.mod["mod1"].obs, "obs key was not copied in mod1"
    assert (
        "Obs_copy" not in output_h5mu.mod["mod2"].obs
    ), "obs key should not have been copied in mod2"
    assert (
        "Obs copy" not in input_h5mu.mod["mod1"].obs
    ), "obs key should not have been copied in input file"
    assert np.all(
        output_h5mu.mod["mod1"].obs["Obs"] == output_h5mu.mod["mod1"].obs["Obs_copy"]
    ), "copied obs column should be identical to original obs column"


def test_copy_index(run_component, random_h5mu_path, input_h5mu, input_h5mu_path):
    output_h5mu_path = random_h5mu_path()

    args = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--modality",
        "mod1",
        "--output_obs_key",
        "Obs_copy",
    ]

    run_component(args)

    assert output_h5mu_path.is_file()

    output_h5mu = mu.read_h5mu(output_h5mu_path)

    assert "Obs_copy" in output_h5mu.mod["mod1"].obs, "obs key was not copied in mod1"
    assert (
        "Obs_copy" not in output_h5mu.mod["mod2"].obs
    ), "obs key should not have been copied in mod2"
    assert (
        "Obs copy" not in input_h5mu.mod["mod1"].obs
    ), "obs key should not have been copied in input file"
    assert np.all(
        output_h5mu.mod["mod1"].obs.index == output_h5mu.mod["mod1"].obs["Obs_copy"]
    ), "copied obs column should be identical to original obs index"


def test_overwrite_keys(run_component, random_h5mu_path, input_h5mu, input_h5mu_path):
    output_h5mu_path = random_h5mu_path()

    args = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--modality",
        "mod1",
        "--input_obs_key",
        "Obs",
        "--output_obs_key",
        "Obs",
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: --output_obs_key already exists: \`Obs\`. Data can not be duplicated.",
        err.value.stdout.decode("utf-8"),
    )

    disable_raise_args = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--modality",
        "mod1",
        "--overwrite_existing_key",
        "--input_obs_key",
        "Obs",
        "--output_obs_key",
        "Obs",
    ]

    run_component(disable_raise_args)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
