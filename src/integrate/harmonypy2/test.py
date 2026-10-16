import sys
import subprocess
import pytest
import mudata
import numpy as np
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "executable": "./target/executable/integrate/harmonypy2/harmonypy2",
    "resources_dir": "./resources_test/pbmc_1k_protein_v3/",
    "config": "src/integrate/harmonypy2/config.vsh.yaml",
}
## VIASH END

input_file = f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def input_with_covariates(write_mudata_to_file):
    input_mudata = mudata.read_h5mu(input_file)
    rna = input_mudata.mod["rna"]
    rng = np.random.default_rng(0)
    rna.obs["batch"] = rng.choice(["a", "b", "c"], size=rna.n_obs)
    rna.obs["lab"] = rng.choice(["x", "y"], size=rna.n_obs)
    return write_mudata_to_file(input_mudata)


@pytest.fixture
def input_with_missing_covariate(write_mudata_to_file):
    input_mudata = mudata.read_h5mu(input_file)
    rna = input_mudata.mod["rna"]
    batch = np.array(
        ["a", "b"] * (rna.n_obs // 2) + ["a"] * (rna.n_obs % 2), dtype=object
    )
    batch[0] = None
    rna.obs["batch"] = batch
    return write_mudata_to_file(input_mudata)


def test_harmonypy2(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_file,
            "--modality",
            "rna",
            "--obsm_input",
            "X_pca",
            "--obsm_output",
            "X_pca_int",
            "--obs_covariates",
            "harmony_integration_leiden_1.0",
            "--output",
            output_path,
            "--output_compression",
            "gzip",
        ]
    )
    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    assert "X_pca_int" in output_data.mod["rna"].obsm
    integrated = output_data.mod["rna"].obsm["X_pca_int"]
    assert integrated.shape == output_data.mod["rna"].obsm["X_pca"].shape
    assert np.all(np.isfinite(integrated))
    assert not np.allclose(integrated, output_data.mod["rna"].obsm["X_pca"])
    # Delete output slots in order to tests for equality with
    # input data
    del output_data["rna"].obsm["X_pca_int"]
    assert_annotation_objects_equal(input_file, output_data)


def test_harmonypy2_multiple_covariates_and_parameters(
    run_component, random_h5mu_path, input_with_covariates
):
    output_path = random_h5mu_path()

    run_component(
        [
            "--input",
            input_with_covariates,
            "--obs_covariates",
            "batch;lab",
            "--theta",
            "2;1",
            "--lambda",
            "1;0.5",
            "--sigma",
            "0.2",
            "--nclust",
            "10",
            "--max_iter_harmony",
            "5",
            "--max_iter_kmeans",
            "10",
            "--epsilon_cluster",
            "0.00001",
            "--epsilon_harmony",
            "0.0001",
            "--batch_prop_cutoff",
            "0.001",
            "--random_state",
            "1",
            "--output",
            output_path,
        ]
    )
    assert output_path.is_file()

    output_data = mudata.read_h5mu(output_path)
    integrated = output_data.mod["rna"].obsm["X_pca_integrated"]
    assert integrated.shape == output_data.mod["rna"].obsm["X_pca"].shape
    assert np.all(np.isfinite(integrated))


def test_harmonypy2_missing_covariate_values_raises(
    run_component, random_h5mu_path, input_with_missing_covariate
):
    output_path = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_with_missing_covariate,
                "--obs_covariates",
                "batch",
                "--output",
                output_path,
            ]
        )
    assert (
        "Covariate column 'batch' contains 1 missing value(s)."
        in err.value.stdout.decode("utf-8")
    )


def test_harmonypy2_unknown_covariate_raises(run_component, random_h5mu_path):
    output_path = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input_file,
                "--obs_covariates",
                "does_not_exist",
                "--output",
                output_path,
            ]
        )
    assert (
        "Covariate column(s) does_not_exist not found in .obs of modality 'rna'."
        in err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
