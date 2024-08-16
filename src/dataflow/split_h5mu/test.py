import sys
from textwrap import dedent
import pytest
import mudata as mu
import anndata as ad
import numpy as np
import pandas as pd
import subprocess
import re


@pytest.fixture
def input_modality_1():
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    obs = pd.DataFrame({'Obs': ["A", "B"], "Obs_nan": [np.nan, np.nan]}, index=df.index)
    var = pd.DataFrame([["a"], ["b"], ["c"]], index=df.columns, columns=["Feat"])
    ad1 = ad.AnnData(df, obs=obs, var=var)
    return ad1


@pytest.fixture
def input_modality_2():
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = ad.AnnData(df, obs=obs2, var=var2)
    return ad2


@pytest.fixture
def input_modality_3():
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    var3 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs3 = pd.DataFrame(["C C", "C_C"], index=df.index, columns=["Obs"])
    ad3 = ad.AnnData(df, obs=obs3, var=var3)
    return ad3


@pytest.fixture
def input_h5mu(input_modality_1, input_modality_2):
    tmp_mudata = mu.MuData({'mod1': input_modality_1, 'mod2': input_modality_2})
    return tmp_mudata


@pytest.fixture
def input_h5mu_path(write_mudata_to_file, input_h5mu):
    return write_mudata_to_file(input_h5mu)


@pytest.fixture
def input_h5mu_non_unique_filenames(input_modality_3):
    tmp_mudata = mu.MuData({'mod3': input_modality_3})
    return tmp_mudata


@pytest.fixture
def input_h5mu_path_non_unique_filenames(write_mudata_to_file, input_h5mu_non_unique_filenames):
    return write_mudata_to_file(input_h5mu_non_unique_filenames)


def test_sample_split(run_component, random_path, input_h5mu, input_h5mu_path):
    output_dir = random_path()
    output_files = random_path(extension="csv")
    args = [
        "--input", input_h5mu_path,
        "--output", str(output_dir),
        "--modality", "mod1",
        "--obs_feature", "Obs",
        "--output_files", str(output_files),
    ]

    run_component(args)
    assert output_files.is_file()
    assert output_dir.is_dir()

    # check output dir and file names
    dir_content = [h5mu_file for h5mu_file in output_dir.iterdir() 
                   if h5mu_file.suffix == ".h5mu" and h5mu_file != input_h5mu_path]
    s1_file = output_dir / f"{input_h5mu_path.stem}_A.h5mu"
    s2_file = output_dir / f"{input_h5mu_path.stem}_B.h5mu"
    assert set(dir_content) == set([s1_file, s2_file])

    # check that number of modalities, variables and observations
    s1 = mu.read_h5mu(s1_file)
    s2 = mu.read_h5mu(s2_file)
    assert s1.n_mod == 2
    assert s2.n_mod == 2

    assert s1.n_obs == input_h5mu.n_obs, "number of observations of split file does not match input file"
    assert s2.n_obs == input_h5mu.n_obs, "number of observations of split file does not match input file"

    assert s1.mod["mod1"].n_obs == 1, "number of observations of split file s1 modality mod1 should equal 1"
    assert s1.mod["mod2"].n_obs == input_h5mu.n_obs, "number of observations of split file s1 modality mod2 should equal input file"

    assert len(s1.mod["mod1"].obs.keys()) == 2, "number of observation keys split file s1 modality mod1 should equal 2"
    assert len(s1.mod["mod2"].obs.keys()) == 1, "number of observation keys split file s1 modality mod2 should equal 1"

    assert s2.mod["mod1"].n_obs == 1, "number of observations of split file s2 modality mod1 should equal 1"
    assert s2.mod["mod2"].n_obs == input_h5mu.n_obs, "number of observations of split file s2 modality mod2 should equal input file"

    assert s1.n_vars == input_h5mu.n_vars, "number of variables of split file s1 should equal input file"
    assert s2.n_vars == input_h5mu.n_vars, "number of variables of split file s1 should equal input file"

    assert s1.mod["mod1"].n_vars == input_h5mu.mod["mod1"].n_vars, "number of variables of split file s1 modalitty mod1 should equal input file"
    assert s1.mod["mod2"].n_vars == input_h5mu.mod["mod1"].n_vars,  "number of variables of split file s1 modalitty mod2 should equal input file"

    assert s2.mod["mod1"].n_vars == input_h5mu.mod["mod1"].n_vars, "number of variables of split file s2 modalitty mod1 should equal input file"
    assert s2.mod["mod2"].n_vars == input_h5mu.mod["mod1"].n_vars,  "number of variables of split file s2 modalitty mod2 should equal input file"

    # check correct sample splitting
    assert np.all(s1.mod["mod1"].obs["Obs"] == "A"), "observation of .obs Obs in s1 should equal A"
    assert np.all(s2.mod["mod1"].obs["Obs"] == "B"), "observation of .obs Obs in s2 should equal B"

    # Check contents of csv file
    expected_csv_output = dedent(
        f"""\
        name,filename
        A,{s1_file.name}
        B,{s2_file.name}
        """
    )
    with open(output_files, 'r') as open_csv_file:
        result = open_csv_file.read()
        assert result == expected_csv_output


def test_sample_split_dropna(run_component, random_path, input_h5mu, input_h5mu_path):
    output_dir = random_path()
    output_files = random_path(extension="csv")
    args = [
        "--input", input_h5mu_path,
        "--output", str(output_dir),
        "--modality", "mod1",
        "--obs_feature", "Obs",
        "--drop_obs_nan", "true",
        "--output_files", str(output_files),
    ]

    run_component(args)

    # check output dir and file names
    s1_file = output_dir / f"{input_h5mu_path.stem}_A.h5mu"
    s2_file = output_dir / f"{input_h5mu_path.stem}_B.h5mu"

    # check that .obs columns with only nan values are dropped correctly
    s1 = mu.read_h5mu(s1_file)
    s2 = mu.read_h5mu(s2_file)

    assert s1.n_obs == input_h5mu.n_obs, "number of observations of split file does not match input file"
    assert s2.n_obs == input_h5mu.n_obs, "number of observations of split file does not match input file"

    assert s1.mod["mod1"].n_obs == 1, "number of observations of split file s1 modality mod1 should equal 1"
    assert s1.mod["mod2"].n_obs == input_h5mu.n_obs, "number of observations of split file s1 modality mod2 should equal input file"

    assert len(s1.mod["mod1"].obs.keys()) == 1, "number of observation keys split file s1 modality mod1 should equal 1"
    assert len(s1.mod["mod2"].obs.keys()) == 1, "number of observation keys split file s1 modality mod2 should equal 1"

def test_sanitizing(run_component, random_path, input_h5mu_path_non_unique_filenames):
    output_dir = random_path()
    output_files = random_path(extension="csv")

    args = [
        "--input", input_h5mu_path_non_unique_filenames,
        "--output", str(output_dir),
        "--modality", "mod3",
        "--obs_feature", "Obs",
        "--drop_obs_nan", "true",
        "--output_files", str(output_files)
    ]

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"ValueError: File names are not unique after sanitizing the --obs_feature Obs values",
        err.value.stdout.decode('utf-8'))

    args_san = [
            "--input", input_h5mu_path_non_unique_filenames,
            "--output", str(output_dir),
            "--modality", "mod3",
            "--obs_feature", "Obs",
            "--drop_obs_nan", "true",
            "--output_files", str(output_files),
            "--ensure_unique_filenames", "true"
        ]

    run_component(args_san)

    # check output dir and file names
    dir_content = [h5mu_file for h5mu_file in output_dir.iterdir()
                   if h5mu_file.suffix == ".h5mu" and h5mu_file != input_h5mu_path_non_unique_filenames]
    s1_file = output_dir / f"{input_h5mu_path_non_unique_filenames.stem}_C_C_0.h5mu"
    s2_file = output_dir / f"{input_h5mu_path_non_unique_filenames.stem}_C_C_1.h5mu"

    assert s1_file.is_file(), f"{s1_file} does not exist"
    assert s2_file.is_file(), f"{s2_file} does not exist"
    assert set(dir_content) == set([s1_file, s2_file]), "Output files do not match file names in csv"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
