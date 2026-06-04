import sys
import pytest
import subprocess
import mudata as mu
import numpy as np
import re

## VIASH START
meta = {
    "name": "foo",
    "resources_dir": "resources_test/",
    "executable": "./target/executable/dimred/pca/pca",
    "config": "./src/dimred/pca/config.vsh.yaml",
}
## VIASH END

input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


@pytest.mark.parametrize("chunked,chunk_size", ([True, 26], [None, None], [None, 26]))
def test_pca(run_component, random_h5mu_path, chunked, chunk_size):
    output_path = random_h5mu_path()
    args = [
        "--input",
        input_path,
        "--output",
        output_path,
        "--obsm_output",
        "X_foo",
        "--num_components",
        "26",
        "--overwrite",
    ]
    if chunked:
        args.append("--chunked")
    if chunk_size:
        args.extend(["--chunk_size", str(chunk_size)])

    # run component
    run_component(args)
    assert output_path.is_file()
    data = mu.read_h5mu(output_path)

    # check whether pca was found
    assert data.mod["rna"].obsm["X_foo"].shape == (data.n_obs, 26)
    assert "highly_variable" not in data.mod["rna"].var.columns
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "pca_variance" in data.mod["rna"].uns
    assert "pca_loadings" in data.mod["rna"].varm
    assert "X_foo" in data.mod["rna"].obsm
    # GH298
    assert not np.array_equal(
        data.mod["rna"].uns["pca_variance"]["variance"],
        data.mod["rna"].uns["pca_variance"]["variance_ratio"],
    )


def test_chunk_size_smaller_than_n_components_raises(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input",
        input_path,
        "--output",
        output_path,
        "--obsm_output",
        "X_foo",
        "--num_components",
        "26",
        "--chunked",
        "--chunk_size",
        "25",
        "--overwrite",
    ]
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"The requested chunk size \(25\) must not be smaller than the number of components \(26\)",
        err.value.stdout.decode("utf-8"),
    )


def test_chunk_size_not_set_raises(run_component, random_h5mu_path):
    output_path = random_h5mu_path()
    args = [
        "--input",
        input_path,
        "--output",
        output_path,
        "--obsm_output",
        "X_foo",
        "--num_components",
        "26",
        "--chunked",
        "--overwrite",
    ]
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(args)
    assert re.search(
        r"Requested to perform an incremental PCA \('chunked'\), but the chunk size is not set.",
        err.value.stdout.decode("utf-8"),
    )


def test_no_overwrite_but_field_also_not_present(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    # create input data
    input_data = mu.read_h5mu(input_path)
    input_data.mod["rna"].uns.pop("pca_variance")
    input_data.mod["rna"].varm.pop("pca_loadings")
    input_data.mod["rna"].obsm.pop("X_pca")
    tmp_file = random_h5mu_path()
    input_data.write(tmp_file)

    # run component
    run_component(
        [
            "--input",
            tmp_file,
            "--output",
            output_path,
            "--obsm_output",
            "X_foo",
            "--num_components",
            "26",
            "--output_compression",
            "gzip",
        ]
    )
    assert output_path.is_file()
    data = mu.read_h5mu(output_path)

    # check whether pca was found
    assert data.mod["rna"].obsm["X_foo"].shape == (data.n_obs, 26)
    assert "highly_variable" not in data.mod["rna"].var.columns
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "pca_variance" in data.mod["rna"].uns
    assert "pca_loadings" in data.mod["rna"].varm
    assert "X_foo" in data.mod["rna"].obsm


def test_selecting_input_layer(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    # generate input data
    input_data = mu.read_h5mu(input_path)
    input_data.mod["rna"].layers["test"] = input_data.mod["rna"].X
    tmp_file = random_h5mu_path()
    input_data.write_h5mu(tmp_file)

    # run component
    run_component(
        [
            "--input",
            tmp_file,
            "--output",
            output_path,
            "--obsm_output",
            "test_foo",
            "--num_components",
            "26",
            "--layer",
            "test",
            "--overwrite",
        ]
    )
    assert output_path.is_file()

    # check whether pca was found
    data = mu.read_h5mu(output_path)
    assert "test_foo" in data.mod["rna"].obsm
    assert data.mod["rna"].obsm["test_foo"].shape == (data.n_obs, 26)
    assert "highly_variable" not in data.mod["rna"].var.columns
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "pca_variance" in data.mod["rna"].uns
    assert "pca_loadings" in data.mod["rna"].varm
    assert "X_pca" in data.mod["rna"].obsm


def test_highly_variable_column_does_not_exist_raises(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                "output.h5mu",
                "--obsm_output",
                "X_foo",
                "--num_components",
                "26",
                "--var_input",
                "does_not_exist",
            ]
        )
    assert (
        "ValueError: Requested to use .var column does_not_exist as "
        "a selection of genes to run the PCA on, but the column is "
        "not available for modality rna" in err.value.stdout.decode("utf-8")
    )


def test_select_highly_variable_column(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    # create input data
    input_data = mu.read_h5mu(input_path)
    input_data.mod["rna"].var["filter_with_hvg"] = True
    tmp_file = random_h5mu_path()
    input_data.write_h5mu(tmp_file)

    # run component
    run_component(
        [
            "--input",
            tmp_file,
            "--output",
            output_path,
            "--obsm_output",
            "X_foo",
            "--num_components",
            "26",
            "--var_input",
            "filter_with_hvg",
            "--overwrite",
        ]
    )
    assert output_path.is_file()

    # check whether pca was found
    data = mu.read_h5mu(output_path)
    assert data.mod["rna"].obsm["X_foo"].shape == (data.n_obs, 26)
    assert "highly_variable" not in data.mod["rna"].var.columns
    assert "filter_with_hvg" in data.mod["rna"].var.columns
    assert "pca_variance" in data.mod["rna"].uns
    assert "pca_loadings" in data.mod["rna"].varm
    assert "X_pca" in data.mod["rna"].obsm


def test_raise_if_input_layer_is_missing(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                "output.h5mu",
                "--obsm_output",
                "X_foo",
                "--num_components",
                "26",
                "--layer",
                "does_not_exist",
                "--var_input",
                "filter_with_hvg",
            ]
        )
    assert (
        "ValueError: does_not_exist was not found in modality rna."
        in err.value.stdout.decode("utf-8")
    )


def test_output_field_already_present_raises(run_component):
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_path,
                "--output",
                "output.h5mu",
                "--obsm_output",
                "X_foo",
                "--num_components",
                "26",
            ]
        )
    assert (
        "ValueError: Requested to create field pca_loadings in .varm for "
        "modality rna, but field already exists." in err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
