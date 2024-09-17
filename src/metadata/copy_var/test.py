import pandas as pd
import numpy as np
import mudata as mu
import anndata as ad
import sys
import pytest


@pytest.fixture
def input_h5mu():
    # generate data
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]], index=df.columns, columns=["Feat"])
    ad1 = ad.AnnData(df, obs=obs, var=var)
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = ad.AnnData(df, obs=obs2, var=var2)
    tmp_mudata = mu.MuData({'mod1': ad1, 'mod2': ad2})
    return tmp_mudata


@pytest.fixture
def input_h5mu_path(write_mudata_to_file, input_h5mu):
    return write_mudata_to_file(input_h5mu)


def test_copy_var(run_component, random_h5mu_path, input_h5mu, input_h5mu_path):
    output_h5mu_path = random_h5mu_path()

    args = [
        "--input", input_h5mu_path,
        "--output", output_h5mu_path,
        "--modality", "mod1",
        "--input_var_key", "Feat",
        "--output_var_key", "Feat_copy"
    ]

    run_component(args)

    assert output_h5mu_path.is_file()

    output_h5mu = mu.read_h5mu(output_h5mu_path)

    assert "Feat_copy" in output_h5mu.mod["mod1"].var, "var key was not copied in mod1"
    assert "Feat_copy" not in output_h5mu.mod["mod2"].var, "var key should not have been copied in mod2"
    assert "Feat_copy" not in input_h5mu.mod["mod1"].var, "var key should not have been copied in input file"
    assert np.all(output_h5mu.mod["mod1"].var["Feat"] == output_h5mu.mod["mod1"].var["Feat_copy"]), "copied var column should be identical to original var column"

def test_copy_index(run_component, random_h5mu_path, input_h5mu, input_h5mu_path):
    output_h5mu_path = random_h5mu_path()

    args = [
        "--input", input_h5mu_path,
        "--output", output_h5mu_path,
        "--modality", "mod1",
        "--output_var_key", "Index_copy"
    ]

    run_component(args)

    assert output_h5mu_path.is_file()

    output_h5mu = mu.read_h5mu(output_h5mu_path)

    assert "Index_copy" in output_h5mu.mod["mod1"].var, "var index was not copied in mod1"
    assert "Index_copy" not in output_h5mu.mod["mod2"].var, "var index should not have been copied in mod2"
    assert "Index_copy" not in input_h5mu.mod["mod1"].var, "var index should not have been copied in input file"
    assert np.all(output_h5mu.mod["mod1"].var.index == output_h5mu.mod["mod1"].var["Index_copy"]), "copied var index should be identical to original var index"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
