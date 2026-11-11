import pytest
import re
import pandas as pd
import uuid
from anndata import AnnData
from mudata import MuData, read_h5mu
from subprocess import CalledProcessError

## VIASH START
meta = {
    "name": "move_obsm_to_obs",
    "resources_dir": "resources_test/",
    "executable": "target/executable/metadata/move_obsm_to_obs/move_obsm_to_obs",
    "config": "src/metadata/move_obsm_to_obs/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def h5mu():
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    obs = pd.DataFrame(
        [["A", "sample1"], ["B", "sample2"]],
        index=df.index,
        columns=["Obs", "sample_id"],
    )
    var = pd.DataFrame(
        [["a", "sample1"], ["b", "sample2"], ["c", "sample1"]],
        index=df.columns,
        columns=["Feat", "sample_id_var"],
    )
    obsm = {
        "obsm_key": pd.DataFrame(
            [["foo", "bar"], ["lorem", "ipsum"]],
            index=obs.index,
            columns=["obsm_col1", "obsm_col2"],
        )
    }
    ad1 = AnnData(df, obs=obs, var=var, obsm=obsm)
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = AnnData(df, obs=obs2, var=var2)
    return MuData({"mod1": ad1, "mod2": ad2})


@pytest.fixture
def write_temp_h5mu(tmp_path):
    def wrapper(test_h5mu):
        test_h5mu_path = tmp_path / f"{str(uuid.uuid4())}.h5mu"
        test_h5mu.write_h5mu(test_h5mu_path)
        return test_h5mu_path

    return wrapper


@pytest.fixture
def h5mu_with_non_overlapping_observations(h5mu):
    h5mu.mod["mod1"].obsm["obsm_key"].index = pd.Index(["obs1", "doesnotexist"])
    return h5mu


def test_move_obsm_to_obs(run_component, h5mu, write_temp_h5mu, tmp_path):
    output = tmp_path / "output.h5mu"
    run_component(
        [
            "--input",
            write_temp_h5mu(h5mu),
            "--modality",
            "mod1",
            "--obsm_key",
            "obsm_key",
            "--output",
            output,
        ]
    )
    assert output.is_file(), "Some output file must have been created."
    output_data = read_h5mu(output)
    pd.testing.assert_index_equal(
        output_data.mod["mod1"].obs.index, pd.Index(["obs1", "obs2"])
    )
    pd.testing.assert_index_equal(
        output_data.mod["mod1"].obs.columns,
        pd.Index(["Obs", "sample_id", "obsm_key_obsm_col1", "obsm_key_obsm_col2"]),
    )
    assert "obsm_key" not in output_data.mod["mod1"].obsm


def test_move_obsm_to_obs_non_overlapping_obs_fails(
    run_component, write_temp_h5mu, h5mu_with_non_overlapping_observations, tmp_path
):
    output = tmp_path / "output.h5mu"
    # Mudata seems to handle this error, but keep this test in just in case mudata drops the ball.
    with pytest.raises((CalledProcessError, ValueError)) as err:
        run_component(
            [
                "--input",
                write_temp_h5mu(h5mu_with_non_overlapping_observations),
                "--modality",
                "mod1",
                "--obsm_key",
                "obsm_key",
                "--output",
                output,
            ]
        )
    expected_message = r"value.index does not match parent’s obs names"
    if isinstance(err, CalledProcessError):
        assert re.search(expected_message, err.value.stdout.decode("utf-8"))
    else:
        assert re.search(expected_message, str(err))


def test_error_non_existing_modality(run_component, h5mu, write_temp_h5mu, tmp_path):
    output = tmp_path / "output.h5mu"
    with pytest.raises(CalledProcessError) as err:
        run_component(
            [
                "--input",
                write_temp_h5mu(h5mu),
                "--modality",
                "foo",
                "--obsm_key",
                "obsm_key",
                "--output",
                output,
            ]
        )
    assert re.search(
        r"ValueError: Modality foo does not exist\.", err.value.stdout.decode("utf-8")
    )


def test_execute_twice_overwrites(run_component, h5mu, write_temp_h5mu, tmp_path):
    output_run_1 = tmp_path / "output1.h5mu"
    run_component(
        [
            "--input",
            write_temp_h5mu(h5mu),
            "--modality",
            "mod1",
            "--obsm_key",
            "obsm_key",
            "--output",
            output_run_1,
        ]
    )
    output_data_run_1 = read_h5mu(output_run_1)
    output_data_run_1.mod["mod1"].obsm = {
        "obsm_key": pd.DataFrame(
            [["dolor", "amet"], ["jommeke", "filiberke"]],
            index=output_data_run_1.mod["mod1"].obs_names,
            columns=["obsm_col1", "obsm_col2"],
        )
    }

    output_run_2 = tmp_path / "output2.h5mu"
    input_run_2 = write_temp_h5mu(output_data_run_1)
    run_component(
        [
            "--input",
            input_run_2,
            "--modality",
            "mod1",
            "--obsm_key",
            "obsm_key",
            "--output",
            output_run_2,
        ]
    )
    assert output_run_2.is_file(), "Some output file must have been created."
    output_data = read_h5mu(output_run_2)
    pd.testing.assert_index_equal(
        output_data.mod["mod1"].obs.index, pd.Index(["obs1", "obs2"])
    )
    pd.testing.assert_index_equal(
        output_data.mod["mod1"].obs.columns,
        pd.Index(["Obs", "sample_id", "obsm_key_obsm_col1", "obsm_key_obsm_col2"]),
    )
    assert "obsm_key" not in output_data.mod["mod1"].obsm


if __name__ == "__main__":
    exit(pytest.main([__file__]))
