import sys
import pytest
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu

## VIASH START
meta = {
    'executable': './target/docker/metadata/add_id/add_id',
    'resources_dir': './src/base/test_utils',
    'cpus': 2,
    'config': './src/metadata/add_id/config.vsh.yaml'
}
## VIASH END


sys.path.append(meta["resources_dir"])
from test_utils.asserters import assert_annotation_objects_equal


def remove_obs_column(mudata_object, column):
    for _, modality in mudata_object.mod.items():
        modality.obs = modality.obs.drop([column],
                                         axis="columns", inplace=False)
    mod_names_columns = [f"{mod_name}:{column}" for mod_name in mudata_object.mod.keys()]
    mudata_object.obs = mudata_object.obs.drop([column] + mod_names_columns,
                                               axis="columns",
                                               inplace=False)
    mudata_object.update()
    return mudata_object

@pytest.fixture
def generate_h5mu():
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
    return tmp_mudata


@pytest.mark.parametrize("output_compression", ["gzip", "lzf", None])
def test_add_id(run_component, generate_h5mu, write_mudata_to_file, random_h5mu_path, output_compression):
    input_data = generate_h5mu
    input_path = write_mudata_to_file(input_data)
    output_path = random_h5mu_path()

    args = [
        "--input", str(input_path),
        "--output", str(output_path),
        "--input_id", "test_id",
    ]
    if output_compression:
        args.extend(["--output_compression", output_compression])
    # run component
    run_component(args)
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "sample_id" in output_data.obs.columns.to_list()
    assert {"test_id"} == set(output_data.obs["sample_id"].to_list())

    # Make sure that the rest of the mudata object stayed the same
    assert_annotation_objects_equal(input_data,
                                    remove_obs_column(output_data, "sample_id"),
                                    check_data=True)


def test_add_id_obs_output(run_component, generate_h5mu,
                           write_mudata_to_file, random_h5mu_path):
    input_data = generate_h5mu
    input_path = write_mudata_to_file(input_data)
    output_path = random_h5mu_path()


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

    # Make sure that the rest of the mudata object stayed the same
    assert_annotation_objects_equal(input_data,
                                    remove_obs_column(output_data, "test_key"),
                                    check_data=True)


def test_add_id_observations_unique(run_component, generate_h5mu,
                                    write_mudata_to_file, random_h5mu_path):
    input_data = generate_h5mu
    input_path = write_mudata_to_file(input_data)
    output_path = random_h5mu_path()

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

    # Make sure that the rest of the mudata object stayed the same
    equals_test_output_data = remove_obs_column(output_data, "sample_id")
    for mod_name, modality in equals_test_output_data.mod.items():
        modality.obs.index =  input_data[mod_name].obs.index
    equals_test_output_data.obs.index = input_data.obs.index
    assert_annotation_objects_equal(input_data, equals_test_output_data, check_data=True)

def test_add_id_overwrites_output_column(run_component, generate_h5mu,
                                         write_mudata_to_file, random_h5mu_path):
    input_data = generate_h5mu
    input_data.obs["already_exists"] = "alread_exists"
    for _, modality in input_data.mod.items():
        modality.obs["already_exists"] = "alread_exists"
    input_path = write_mudata_to_file(input_data)
    output_path = random_h5mu_path()

    # run component
    run_component([
        "--input", str(input_path),
        "--output", str(output_path),
        "--input_id", "test_id",
        "--obs_output", "already_exists"
    ])
    output_data = read_h5mu(output_path)
    assert "already_exists" in output_data.obs.columns.to_list()
    for mod_data in output_data.mod.values():
        assert "already_exists" in mod_data.obs.columns.to_list()
    assert {"test_id"} == set(output_data.obs["already_exists"].to_list())


if __name__ == "__main__":
    sys.exit(pytest.main([__file__], plugins=["test_utils.fixtures"]))


