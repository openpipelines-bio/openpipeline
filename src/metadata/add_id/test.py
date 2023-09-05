import sys
import pytest
import pandas as pd
from anndata import AnnData
from tempfile import NamedTemporaryFile
from mudata import MuData, read_h5mu
from pathlib import Path

## VIASH START
meta = {
    'executable': './target/docker/metadata/add_id/add_id',
    'resources_dir': './resources_test/concat_test_data/',
    'cpus': 2
}
## VIASH END

@pytest.fixture
def generate_h5mu(tmp_path):
    tmp_file = tmp_path / "test.h5mu"

    # generate data
    df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"])
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]],
                        index=df.columns, columns=["Feat"])
    ad1 = AnnData(df, obs=obs, var=var)
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = AnnData(df, obs=obs2, var=var2)
    tmp_mudata = MuData({'mod1': ad1, 'mod2': ad2})
    
    # write to file
    tmp_mudata.write_h5mu(tmp_file)
    return tmp_file, tmp_mudata

def test_add_id(run_component, generate_h5mu, tmp_path):
    input_path, input_data = generate_h5mu
    output_path = tmp_path / "with_id.h5mu"

    # run component
    run_component([
        "--input", str(input_path),
        "--output", str(output_path),
        "--input_id", "test_id",
        "--output_compression", "gzip"
    ])
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "sample_id" in output_data.obs.columns.to_list()
    assert {"test_id"} == set(output_data.obs["sample_id"].to_list())
    pd.testing.assert_index_equal(output_data.obs.index, input_data.obs.index)

def test_add_id_obs_output(run_component, generate_h5mu, tmp_path):
    input_path, input_data = generate_h5mu
    output_path = tmp_path / "with_id.h5mu"

    # run component
    run_component([
        "--input", str(input_path),
        "--output", str(output_path),
        "--input_id", "test_id",
        "--obs_output", "test_key"
    ])
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "test_key" in output_data.obs.columns.to_list()
    for mod_data in output_data.mod.values():
        assert "test_key" in mod_data.obs.columns.to_list()
    assert {"test_id"} == set(output_data.obs["test_key"].to_list())
    pd.testing.assert_index_equal(output_data.obs.index, input_data.obs.index)


def test_add_id_observations_unique(run_component, generate_h5mu, tmp_path):
    input_path, input_data = generate_h5mu
    output_path = tmp_path / "with_id.h5mu"

    # run component
    run_component([
        "--input", str(input_path),
        "--output", str(output_path),
        "--input_id", "test_id",
        "--make_observation_keys_unique"
    ])
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "sample_id" in output_data.obs.columns.to_list()
    assert {"test_id"} == set(output_data.obs["sample_id"].to_list())
    pd.testing.assert_index_equal(output_data.obs.index, "test_id_" + input_data.obs.index)

if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))

