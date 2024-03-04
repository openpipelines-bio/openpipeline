import sys
import pytest
import pandas as pd
from anndata import AnnData
from mudata import MuData, read_h5mu
from openpipelinetestutils.asserters import assert_annotation_objects_equal
from openpipelinetestutils.utils import remove_annotation_column

## VIASH START
meta = {
    'executable': './target/docker/metadata/add_id/add_id',
    'resources_dir': './resources_test/concat_test_data/',
    'cpus': 2,
    'config': './src/metadata/add_id/config.vsh.yaml'
}
## VIASH END

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
def test_add_id(run_component, small_mudata, small_mudata_path, random_h5mu_path, output_compression):
    output_path = random_h5mu_path()

    args = [
        "--input", str(small_mudata_path),
        "--output", str(output_path),
        "--input_id", "test_id",
    ]

# TODO order of columns not important
    expected_obs_mod1 = pd.DataFrame(
        {
            "Obs": ["A", "B"],
            "sample_id": ["test_id", "test_id"]
        },
        index = pd.Index(["obs1", "obs2"])
    ).astype(
        {
            "Obs": "object",
            "sample_id": "category"
        }
    )
    expected_obs_mod2 = pd.DataFrame(
        {
            "Obs": ["C", "D"],
            "sample_id": ["test_id", "test_id"]
        },
        index = pd.Index(["obs1", "obs2"])
    ).astype(
        {
            "Obs": "object",
            "sample_id": "category"
        }
    )
    expected_obs_mudata = pd.DataFrame(
        {
            "mod1:Obs": ["A", "B"],
            "mod1:sample_id": ["test_id", "test_id"],
            "mod2:Obs": ["C", "D"],
            "mod2:sample_id": ["test_id", "test_id"],
            "sample_id": ["test_id", "test_id"]
        },
        index = pd.Index(["obs1", "obs2"])
    ).astype(
        {
            "mod1:Obs": "object",
            "mod1:sample_id": "category",
            "mod2:Obs": "object",
            "mod2:sample_id": "category",
            "sample_id": "category"
        }
    )
    small_mudata.mod["mod1"].obs = expected_obs_mod1
    small_mudata.mod["mod2"].obs = expected_obs_mod2
    small_mudata.obs = expected_obs_mudata

    if output_compression:
        args.extend(["--output_compression", output_compression])
    # run component
    run_component(args)
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)

    assert_annotation_objects_equal(output_data, small_mudata)


def test_add_id_obs_output(run_component, generate_h5mu,
                           write_mudata_to_file, random_h5mu_path):
    small_mudata = generate_h5mu
    small_mudata_path = write_mudata_to_file(small_mudata)
    output_path = random_h5mu_path()


    # run component
    run_component([
        "--input", str(small_mudata_path),
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
    pd.testing.assert_index_equal(output_data.obs.index, small_mudata.obs.index)

    # Make sure that the rest of the mudata object stayed the same
    assert_annotation_objects_equal(small_mudata,
                                    remove_annotation_column(output_data, "test_key", "obs"),
                                    check_data=True)


def test_add_id_observations_unique(run_component, generate_h5mu,
                                    write_mudata_to_file, random_h5mu_path):
    small_mudata = generate_h5mu
    small_mudata_path = write_mudata_to_file(small_mudata)
    output_path = random_h5mu_path()

    # run component
    run_component([
        "--input", str(small_mudata_path),
        "--output", str(output_path),
        "--input_id", "test_id",
        "--make_observation_keys_unique"
    ])
    assert output_path.is_file()

    # check output
    output_data = read_h5mu(output_path)
    assert "sample_id" in output_data.obs.columns.to_list()
    assert {"test_id"} == set(output_data.obs["sample_id"].to_list())
    pd.testing.assert_index_equal(output_data.obs.index, "test_id_" + small_mudata.obs.index)

    # Make sure that the rest of the mudata object stayed the same
    equals_test_output_data = remove_annotation_column(output_data, "sample_id", "obs")
    for mod_name, modality in equals_test_output_data.mod.items():
        modality.obs.index =  small_mudata[mod_name].obs.index
    equals_test_output_data.obs.index = small_mudata.obs.index
    assert_annotation_objects_equal(small_mudata, equals_test_output_data, check_data=True)

def test_add_id_overwrites_output_column(run_component, generate_h5mu,
                                         write_mudata_to_file, random_h5mu_path):
    small_mudata = generate_h5mu
    small_mudata.obs["already_exists"] = "alread_exists"
    for _, modality in small_mudata.mod.items():
        modality.obs["already_exists"] = "alread_exists"
    small_mudata_path = write_mudata_to_file(small_mudata)
    output_path = random_h5mu_path()

    # run component
    run_component([
        "--input", str(small_mudata_path),
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
    sys.exit(pytest.main([__file__]))


