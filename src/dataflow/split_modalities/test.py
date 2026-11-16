import sys
import pytest
import mudata as md
import anndata as ad
import pandas as pd
from io import StringIO
from openpipeline_testutils.asserters import assert_annotation_objects_equal
from textwrap import dedent

## VIASH START
meta = {
    "name": "./target/native/dataflow/split_modalities/split_modalities",
    "resources_dir": "./resources_test/",
    "config": "./src/dataflow/split_modalities/config.vsh.yaml",
    "executable": "./target/docker/dataflow/split_modalities/split_modalities",
}
## VIASH END


@pytest.fixture(params=["h5", "zarr"])
def file_format(request):
    return request.param


@pytest.fixture
def input_modality_1():
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    obs = pd.DataFrame([["A"], ["B"]], index=df.index, columns=["Obs"])
    var = pd.DataFrame([["a"], ["b"], ["c"]], index=df.columns, columns=["Feat"])
    ad1 = ad.AnnData(df, obs=obs, var=var)
    return ad1


@pytest.fixture
def input_modality_2():
    df = pd.DataFrame(
        [[1, 2, 3], [4, 5, 6]], index=["obs1", "obs2"], columns=["var1", "var2", "var3"]
    )
    var2 = pd.DataFrame(["d", "e", "g"], index=df.columns, columns=["Feat"])
    obs2 = pd.DataFrame(["C", "D"], index=df.index, columns=["Obs"])
    ad2 = ad.AnnData(df, obs=obs2, var=var2)
    return ad2


@pytest.fixture
def input_h5mu(input_modality_1, input_modality_2):
    tmp_mudata = md.MuData({"mod1": input_modality_1, "mod2": input_modality_2})
    return tmp_mudata


@pytest.fixture
def input_h5mu_path(input_h5mu, random_path, file_format):
    output_path = random_path("h5mu" if file_format != "zarr" else "zarr")
    write_func = input_h5mu.write_zarr if file_format == "zarr" else input_h5mu.write
    write_func(output_path)
    return output_path


@pytest.mark.parametrize("compression", ["gzip", None])
def test_split(
    run_component,
    random_path,
    input_h5mu,
    input_h5mu_path,
    input_modality_1,
    input_modality_2,
    compression,
    file_format,
):
    output_dir = random_path()
    output_types = random_path(extension="csv")
    args = [
        "--input",
        input_h5mu_path,
        "--output",
        str(output_dir),
        "--output_types",
        str(output_types),
    ]
    if compression:
        args += ["--output_compression", compression]
    run_component(args)
    assert output_types.is_file()
    assert output_dir.is_dir()
    output_extension = "h5mu" if file_format != "zarr" else "zarr"
    # check output dir
    dir_content = [
        h5mu_file
        for h5mu_file in output_dir.iterdir()
        if h5mu_file.suffix == f".{output_extension}" and h5mu_file != input_h5mu_path
    ]
    mod1_file = output_dir / f"{input_h5mu_path.stem}_mod1.{output_extension}"
    mod2_file = output_dir / f"{input_h5mu_path.stem}_mod2.{output_extension}"
    assert set(dir_content) == set([mod1_file, mod2_file])
    reader = md.read_zarr if file_format == "zarr" else md.read_h5mu
    mod1 = reader(mod1_file)
    mod2 = reader(mod2_file)
    assert mod1.n_mod == 1
    assert mod2.n_mod == 1

    assert_annotation_objects_equal(mod1.mod["mod1"], input_modality_1)
    assert_annotation_objects_equal(mod2.mod["mod2"], input_modality_2)

    assert mod1.n_obs == input_h5mu.n_obs
    assert mod2.n_obs == input_h5mu.n_obs

    assert mod1.var_names.equals(input_h5mu.mod["mod1"].var_names)
    assert mod2.var_names.equals(input_h5mu.mod["mod2"].var_names)

    # mudata >= 0.4 does not pull modality columns into the global .var,
    # so the annotation is only checked at the modality level.
    assert set(mod1.mod["mod1"].var.columns) == set(input_h5mu.mod["mod1"].var.columns)
    assert set(mod2.mod["mod2"].var.columns) == set(input_h5mu.mod["mod2"].var.columns)

    expected_csv_output = dedent(
        f"""\
        name,filename
        mod1,{mod1_file.name}
        mod2,{mod2_file.name}
        """
    )
    expected_csv = pd.read_csv(
        StringIO(expected_csv_output), sep=",", header=0, index_col=0
    )
    found_csv = pd.read_csv(output_types, sep=",", header=0, index_col=0)
    pd.testing.assert_frame_equal(expected_csv, found_csv, check_like=True)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
