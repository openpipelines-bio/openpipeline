import sys
import pytest
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu

## VIASH START
meta = {
    "executable": "./target/executable/transform/move_layer/move_layer",
    "config": "./src/transform/move_layer/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def test_mudata(random_h5mu_path):
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
    obsm = pd.DataFrame(
        [["X", "W"]], index=pd.Index([0]), columns=["uns_col1", "uns_col2"]
    )
    ad1 = AnnData(df, obs=obs, var=var, uns={"obsm1": obsm})
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = AnnData(df, obs=obs2, var=var2)

    test_h5mu = random_h5mu_path()
    mudata = MuData({"mod1": ad1, "mod2": ad2})
    mudata.write_h5mu(test_h5mu)
    return test_h5mu


def test_move_layer(test_mudata, run_component, random_h5mu_path):
    output_file = random_h5mu_path()
    run_component(
        [
            "--input",
            test_mudata,
            "--modality",
            "mod1",
            "--output_layer",
            "test_layer",
            "--output",
            output_file,
        ]
    )
    assert output_file.is_file()
    output_mudata = read_h5mu(output_file)
    assert "test_layer" in output_mudata.mod["mod1"].layers
    assert output_mudata.mod["mod1"].X is None


def test_move_layer_select_input_layer(test_mudata, run_component, random_h5mu_path):
    output_file = random_h5mu_path()
    run_component(
        [
            "--input",
            test_mudata,
            "--modality",
            "mod1",
            "--output_layer",
            "test_layer",
            "--output",
            output_file,
        ]
    )
    output_file_2 = random_h5mu_path()
    run_component(
        [
            "--input",
            output_file,
            "--modality",
            "mod1",
            "--input_layer",
            "test_layer",
            "--output_layer",
            "test_layer2",
            "--output",
            output_file_2,
        ]
    )
    assert output_file_2.is_file()
    output_mudata = read_h5mu(output_file_2)
    assert "test_layer2" in output_mudata.mod["mod1"].layers
    assert "test_layer" not in output_mudata.mod["mod1"].layers
    assert output_mudata.mod["mod1"].X is None


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
