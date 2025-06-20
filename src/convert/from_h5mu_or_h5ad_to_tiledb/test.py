import sys
import pytest
import re
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import tiledbsoma
import pyarrow as pa
from subprocess import CalledProcessError

ad.settings.allow_write_nullable_strings = True


## VIASH START
meta = {
    "resources_dir": "resources_test",
    "executable": "./target/executable/convert/from_h5mu_or_h5ad_to_tiledb/from_h5mu_or_h5ad_to_tiledb",
    "config": "./src/convert/from_h5mu_or_h5ad_to_tiledb/config.vsh.yaml",
}
## VIASH END


@pytest.fixture
def sample_1_modality_1():
    """
    >>> ad1.obs
          Obs1 Shared_obs_col
     obs1    A              B
     obs2    C              D

    >>> ad1.var
                         Feat1  Shared_feat_col
    var1                     a                b
    var2                     c                d
    overlapping_var_mod1     e                f

    >>> ad1.X
    array([[1, 2, 3],
           [4, 5, 6]])
    """

    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]],
        index=["obs1", "obs2"],
        columns=["var1", "var2", "overlapping_var_mod1"],
    )

    layer_1 = np.array([[94, 81, 60], [50, 88, 69]])

    obs = pd.DataFrame(
        [["A", "B", True, 0.9, 1, pd.NA], ["C", "D", False, 0.2, 2, 8]],
        index=df.index,
        columns=[
            "Obs1",
            "Shared_obs_col",
            "bool_column",
            "float_column",
            "int_column",
            "Shared_obs_col_int",
        ],
    ).astype(
        {
            "Obs1": "object",
            "Shared_obs_col": "object",
            "bool_column": "bool",
            "float_column": "float64",
            "int_column": "int64",
            "Shared_obs_col_int": "Int64",
        }
    )

    var = pd.DataFrame(
        [["a", "b"], ["c", "d"], ["e", "f"]],
        index=df.columns,
        columns=["Feat1", "Shared_feat_col"],
    )
    varm = np.random.rand(df.columns.size, 5)
    varm2 = pd.DataFrame(
        np.random.rand(df.columns.size, 5),
        index=var.index,
        columns=["varmcol1", "varmcol2", "varmcol3", "varmcol4", "varmcol5"],
    )
    ad1 = ad.AnnData(
        X=df,
        layers={"layer1": layer_1},
        obs=obs,
        var=var,
        varm={"random_vals_mod1": varm, "df_varm": varm2},
        uns={
            "uns_unique_to_sample1": pd.DataFrame(
                ["foo"], index=["bar"], columns=["col1"]
            ),
            "overlapping_uns_key": pd.DataFrame(
                ["jing"], index=["jang"], columns=["col2"]
            ),
        },
    )
    return ad1


@pytest.fixture
def sample_1_modality_2():
    """
    >>> ad2.X
    array([[ 7,  8],
           [ 9, 10],
           [11, 12]])

    >>> ad2.obs
         Obs2 Obs3 Shared_obs_col
    obs3    E    F              G
    obs2    H    I              J
    obs5    K    L              M

    >>> ad2.var
                           Feat2  Shared_feat_col
    overlapping_var_mod1       d                e
    var5                       f                g

    """
    df = pd.DataFrame(
        [[7, 8], [9, 10], [11, 12]],
        index=["obs3", "obs2", "obs5"],
        columns=["overlapping_var_mod1", "var4"],
    )

    layer_1 = np.array([[20, 35], [76, 93], [100, 38]])

    obs = pd.DataFrame(
        [["E", "F", "G", 1], ["H", "I", "J", 2.0], ["K", "L", "M", 3.0]],
        index=df.index,
        columns=["Obs2", "Obs3", "Shared_obs_col", "Shared_obs_col_int"],
    )
    var = pd.DataFrame(
        [["d", "e"], ["f", "g"]], index=df.columns, columns=["Feat2", "Shared_feat_col"]
    )
    ad2 = ad.AnnData(df, obs=obs, var=var, layers={"layer_foo": layer_1})
    return ad2


@pytest.fixture
def sample_1_mudata(sample_1_modality_1, sample_1_modality_2):
    return md.MuData({"mod1": sample_1_modality_1, "mod2": sample_1_modality_2})


@pytest.fixture
def random_h5ad_path(random_path):
    def wrapper():
        return random_path(extension="h5ad")

    return wrapper


@pytest.fixture
def sample_1_mudata_file(sample_1_mudata, random_h5mu_path):
    output_path = random_h5mu_path()
    sample_1_mudata.write(output_path)
    return output_path


@pytest.fixture
def sample_1_modality_1_anndata_file(sample_1_modality_1, random_h5ad_path):
    output_path = random_h5ad_path()
    sample_1_modality_1.write(output_path)
    return output_path


def test_convert_anndata(run_component, sample_1_modality_1_anndata_file, tmp_path):
    output_path = tmp_path / "output"
    run_component(
        [
            "--input",
            sample_1_modality_1_anndata_file,
            "--rna_modality",
            "mod1",
            "--rna_raw_layer_input",
            "X",
            "--rna_normalized_layer_input",
            "layer1",
            "--rna_var_gene_names_input",
            "Feat1",
            "--tiledb_dir",
            output_path,
        ]
    )
    with tiledbsoma.DataFrame.open(f"{output_path}/obs") as obs_arr:
        obs_contents = obs_arr.read().concat()
        expected_schema = pa.schema(
            [
                pa.field("soma_joinid", pa.int64(), nullable=False),
                pa.field("cell_id", pa.large_string()),
                pa.field("Obs1", pa.large_string()),
                pa.field("Shared_obs_col", pa.large_string()),
                pa.field("bool_column", pa.bool_()),
                pa.field("float_column", pa.float64()),
                pa.field("int_column", pa.int64()),
                pa.field("Shared_obs_col_int", pa.int64()),
            ]
        )
        expected = pa.Table.from_arrays(
            [
                pa.array([0, 1]),
                pa.array(["obs1", "obs2"]),
                pa.array(["A", "C"]),
                pa.array(["B", "D"]),
                pa.array([True, False]),
                pa.array([0.9, 0.2]),
                pa.array([1, 2]),
                pa.array([None, 8]),
            ],
            schema=expected_schema,
        )
        assert expected.equals(obs_contents)

    with tiledbsoma.DataFrame.open(f"{output_path}/ms/rna/var") as obs_arr:
        var_contents = obs_arr.read().concat()
        expected_schema = pa.schema(
            [
                pa.field("soma_joinid", pa.int64(), nullable=False),
                pa.field("rna_index", pa.large_string()),
                pa.field("Feat1", pa.large_string()),
                pa.field("Shared_feat_col", pa.large_string()),
                pa.field("gene_symbol", pa.large_string()),
            ]
        )
        expected = pa.Table.from_arrays(
            [
                pa.array([0, 1, 2]),
                pa.array(["var1", "var2", "overlapping_var_mod1"]),
                pa.array(["a", "c", "e"]),
                pa.array(["b", "d", "f"]),
                pa.array(["Feat1", "Feat1", "Feat1"]),
            ],
            schema=expected_schema,
        )
        assert expected.equals(var_contents)
        with tiledbsoma.Experiment.open(f"{output_path}") as exp:
            with tiledbsoma.ExperimentAxisQuery(
                experiment=exp, measurement_name="rna"
            ) as query:
                assert query.n_obs == 2
                data = query.to_anndata("raw").X.todense()
                expected = (
                    [
                        [1, 2, 3],
                        [4, 5, 6],
                    ],
                )
                np.testing.assert_allclose(data, np.asmatrix(np.array(expected)))
                log_normalized_data = query.to_anndata("log_normalized").X.todense()
                expected_log_normalized = [
                    [94, 81, 60],
                    [50, 88, 69],
                ]
                np.testing.assert_allclose(
                    log_normalized_data, np.asmatrix(np.array(expected_log_normalized))
                )


@pytest.mark.parametrize(
    "extra_arg",
    [
        ("--prot_modality", "mod2"),
        ("--prot_raw_layer_input", "X"),
        ("--prot_normalized_layer_input", "layer_foo"),
    ],
)
def test_convert_anndata_extra_prot_arguments_raises(
    run_component, sample_1_modality_1_anndata_file, tmp_path, extra_arg
):
    output_path = tmp_path / "output"
    args = [
        "--input",
        sample_1_modality_1_anndata_file,
        "--rna_modality",
        "mod1",
        "--rna_raw_layer_input",
        "X",
        "--rna_normalized_layer_input",
        "layer1",
        "--rna_var_gene_names_input",
        "Feat1",
        "--tiledb_dir",
        output_path,
    ]
    args += extra_arg
    with pytest.raises(CalledProcessError) as err:
        run_component(args)
    assert re.search(
        rf"'{extra_arg[0].lstrip('--')}' can not be used when the input is not a MuData file\.",
        err.value.stdout.decode("utf-8"),
    )


@pytest.mark.parametrize(
    "params",
    [
        [("--rna_raw_layer_input", "X"), ("--rna_normalized_layer_input", "X")],
        [("--rna_raw_layer_output", "X"), ("--rna_normalized_layer_output", "X")],
        [("--rna_modality", "rna"), ("--prot_modality", "rna")],
        [("--prot_raw_layer_input", "X"), ("--prot_normalized_layer_input", "X")],
        [("--prot_raw_layer_output", "X"), ("--prot_normalized_layer_output", "X")],
    ],
)
def test_arg_values_should_not_be_the_same(
    run_component, sample_1_mudata_file, tmp_path, params
):
    output_path = tmp_path / "output"
    args = [
        "--input",
        sample_1_mudata_file,
        "--rna_var_gene_names_input",
        "Feat1",
        "--tiledb_dir",
        output_path,
    ]
    args += [param for param_tup in params for param in param_tup]
    # The following avoids complaints about misszing input argument
    if "--rna_raw_layer_input" not in args:
        args += [
            "--rna_raw_layer_input",
            "X",
            "--rna_normalized_layer_input",
            "log_normalized",
        ]
    if "--rna_modality" not in args:
        args += ["--rna_modality", "mod1"]

    with pytest.raises(CalledProcessError) as err:
        run_component(args)
    assert re.search(
        rf"The value for argument '{params[0][0].lstrip('-')}' \({params[0][1]}\) must not be the same as for argument '{params[1][0].lstrip('-')}' \({params[1][1]}\)",
        err.value.stdout.decode("utf-8"),
    )


def test_output_directory_already_exists(run_component, sample_1_mudata_file, tmp_path):
    output_path = tmp_path / "output"
    output_path.mkdir()
    run_component(
        [
            "--input",
            sample_1_mudata_file,
            "--rna_modality",
            "mod1",
            "--prot_modality",
            "mod2",
            "--rna_raw_layer_input",
            "X",
            "--rna_normalized_layer_input",
            "layer1",
            "--prot_raw_layer_input",
            "X",
            "--prot_normalized_layer_input",
            "layer_foo",
            "--rna_var_gene_names_input",
            "Feat1",
            "--tiledb_dir",
            output_path,
        ]
    )


def test_convert_mudata(
    run_component, sample_1_mudata_file, sample_1_modality_1, tmp_path
):
    output_path = tmp_path / "output"
    run_component(
        [
            "--input",
            sample_1_mudata_file,
            "--rna_modality",
            "mod1",
            "--prot_modality",
            "mod2",
            "--rna_raw_layer_input",
            "X",
            "--rna_normalized_layer_input",
            "layer1",
            "--prot_raw_layer_input",
            "X",
            "--prot_normalized_layer_input",
            "layer_foo",
            "--rna_var_gene_names_input",
            "Feat1",
            "--tiledb_dir",
            output_path,
        ]
    )

    assert output_path.is_dir()
    with tiledbsoma.DataFrame.open(f"{output_path}/obs") as obs_arr:
        obs_contents = obs_arr.read().concat()
        expected_schema = pa.schema(
            [
                pa.field("soma_joinid", pa.int64(), nullable=False),
                pa.field("cell_id", pa.large_string()),
                pa.field("Obs1", pa.large_string()),
                pa.field("Shared_obs_col", pa.large_string()),
                pa.field("bool_column", pa.bool_()),
                pa.field("float_column", pa.float64()),
                pa.field("int_column", pa.int64()),
                pa.field("Shared_obs_col_int", pa.float64()),
                pa.field("Obs2", pa.large_string()),
                pa.field("Obs3", pa.large_string()),
            ]
        )
        expected = pa.Table.from_arrays(
            [
                pa.array([0, 1, 2, 3]),
                pa.array(["obs1", "obs2", "obs3", "obs5"]),
                pa.array(["A", "C", None, None]),
                pa.array(["B", "D", "G", "M"]),
                pa.array([True, False, None, None]),
                pa.array([0.9, 0.2, None, None]),
                pa.array([1, 2, None, None]),
                pa.array([None, 8, 1, 3]),
                pa.array([None, None, "E", "K"]),
                pa.array([None, None, "F", "L"]),
            ],
            schema=expected_schema,
        )
        # Convert to pandas in order to allow the order of the fields to differ
        pd.testing.assert_frame_equal(
            expected.to_pandas(), obs_contents.to_pandas(), check_like=True
        )

    with tiledbsoma.DataFrame.open(f"{output_path}/ms/rna/var") as obs_arr:
        var_contents = obs_arr.read().concat()
        expected_schema = pa.schema(
            [
                pa.field("soma_joinid", pa.int64(), nullable=False),
                pa.field("rna_index", pa.large_string()),
                pa.field("Feat1", pa.large_string()),
                pa.field("Shared_feat_col", pa.large_string()),
                pa.field("gene_symbol", pa.large_string()),
            ]
        )
        expected = pa.Table.from_arrays(
            [
                pa.array([0, 1, 2]),
                pa.array(["var1", "var2", "overlapping_var_mod1"]),
                pa.array(["a", "c", "e"]),
                pa.array(["b", "d", "f"]),
                pa.array(["Feat1", "Feat1", "Feat1"]),
            ],
            schema=expected_schema,
        )
        assert expected.equals(var_contents)
    with tiledbsoma.Experiment.open(f"{output_path}") as exp:
        with tiledbsoma.ExperimentAxisQuery(
            experiment=exp, measurement_name="rna"
        ) as query:
            assert query.n_obs == 4
            data = query.to_anndata("raw").X.todense()
            expected = ([[1, 2, 3], [4, 5, 6], [0, 0, 0], [0, 0, 0]],)
            np.testing.assert_allclose(data, np.asmatrix(np.array(expected)))
            log_normalized_data = query.to_anndata("log_normalized").X.todense()
            expected_log_normalized = [[94, 81, 60], [50, 88, 69], [0, 0, 0], [0, 0, 0]]
            np.testing.assert_allclose(
                log_normalized_data, np.asmatrix(np.array(expected_log_normalized))
            )

        with tiledbsoma.ExperimentAxisQuery(
            experiment=exp, measurement_name="prot"
        ) as query:
            assert query.n_obs == 4
            data = query.to_anndata("raw").X.todense()
            expected = ([[0, 0], [9, 10], [7, 8], [11, 12]],)
            np.testing.assert_allclose(data, np.asmatrix(np.array(expected)))
            log_normalized_data = query.to_anndata("log_normalized").X.todense()
            expected_log_normalized = [[0, 0], [76, 93], [20, 35], [100, 38]]
            np.testing.assert_allclose(
                log_normalized_data, np.asmatrix(np.array(expected_log_normalized))
            )
    with tiledbsoma.Collection.open(
        f"{output_path}/ms/rna/varm", "r"
    ) as multidim_collection:
        assert set(multidim_collection.keys()) == {"df_varm", "random_vals_mod1"}
        for varm_key in multidim_collection:
            expected_data = sample_1_modality_1.varm[varm_key]
            try:
                expected = expected_data.to_numpy()
                expected_index = sample_1_modality_1.varm[varm_key].columns.to_list()
                col_index = [
                    multidim_collection[varm_key].metadata[str(index_)]
                    for index_ in range(len(expected_index))
                ]
                assert expected_index == col_index
            except AttributeError:
                expected = expected_data
            contents = (
                multidim_collection[varm_key]
                .read()
                .coos()
                .concat()
                .to_scipy()
                .todense()
            )
            np.testing.assert_allclose(expected, contents)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
