from pathlib import Path
import re
import subprocess
import pytest

import mudata as mu
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix, csr_array

## VIASH START
meta = {
    "name": "foo",
    "resources_dir": "resources_test/",
    "executable": "target/executable/filter/filter_with_scrublet/filter_with_scrublet",
    "config": "./src/filter/filter_with_scrublet/config.vsh.yaml",
}
## VIASH END

# read input file


@pytest.fixture
def input_mudata_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def input_mudata(input_mudata_path):
    return mu.read_h5mu(input_mudata_path)


@pytest.fixture
def input_with_failed_run(random_h5mu_path, input_mudata_path):
    new_mudata_path = random_h5mu_path()

    mudata_in = mu.read_h5mu(input_mudata_path)

    # Make test reproducable
    np.random.seed(4)

    # Simulate a failed scrublet run by passing very little cells
    mudata = mudata_in[152].copy()
    nobs = 100
    x_data = np.repeat(mudata.mod["rna"].X.todense(), nobs, axis=0)

    # Random perturbations because otherwise the detection fails in other ways (PCA cannot be run)
    replace_rate = 0.000001
    mask = np.random.choice(
        [0, 1], size=x_data.shape, p=((1 - replace_rate), replace_rate)
    ).astype("bool")
    r = np.random.rand(*x_data.shape) * np.max(x_data)
    x_data[mask] = r[mask]

    # create obs
    obs_name = mudata.mod["rna"].obs.index.to_list()[0]
    obs_data = pd.DataFrame([], index=[f"{obs_name}_{i}" for i in range(nobs)])

    # create resulting mudata
    mod = ad.AnnData(X=csr_matrix(x_data), obs=obs_data, var=mudata.mod["rna"].var)
    new_mudata = mu.MuData({"rna": mod})
    new_mudata.update()
    new_mudata.write(new_mudata_path)

    return new_mudata_path


def test_filter_a_little_bit(
    run_component, random_h5mu_path, input_mudata_path, input_mudata
):
    output_mu = random_h5mu_path()

    run_component(
        [
            "--input",
            input_mudata_path,
            "--output",
            output_mu,
            "--min_counts",
            "3",
            "--output_compression",
            "gzip",
        ]
    )
    assert Path(output_mu).is_file(), "Output file not found"

    mu_out = mu.read_h5mu(output_mu)
    assert "filter_with_scrublet" in mu_out.mod["rna"].obs

    new_obs = mu_out.mod["rna"].n_obs
    new_vars = mu_out.mod["rna"].n_vars
    assert (
        new_obs == input_mudata.mod["rna"].n_obs
    ), "No RNA obs should have been filtered"
    assert (
        new_vars == input_mudata.mod["rna"].n_vars
    ), "No RNA vars should have been filtered"
    assert (
        mu_out.mod["prot"].n_obs == input_mudata.mod["prot"].n_obs
    ), "No prot obs should have been filtered"
    assert (
        mu_out.mod["prot"].n_vars == input_mudata.mod["prot"].n_vars
    ), "No prot vars should have been filtered"
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ], "Feature types of RNA modality should be Gene Expression"
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ], "Feature types of prot modality should be Antibody Capture"


def test_filtering_a_lot(
    run_component, random_h5mu_path, input_mudata_path, input_mudata
):
    output_mu = random_h5mu_path()

    run_component(
        [
            "--input",
            input_mudata_path,
            "--output",
            output_mu,
            "--modality",
            "rna",
            "--min_counts",
            "10",
            "--num_pca_components",
            "10",
            "--do_subset",
        ]
    )
    assert Path(output_mu).is_file(), "Output file not found"

    mu_out = mu.read_h5mu(output_mu)
    new_obs = mu_out.mod["rna"].n_obs
    new_vars = mu_out.mod["rna"].n_vars
    assert (
        new_obs < input_mudata.mod["rna"].n_obs
    ), "Some cells should have been filtered"
    assert (
        new_vars == input_mudata.mod["rna"].n_vars
    ), "No genes should have been filtered"
    assert (
        mu_out.mod["prot"].n_obs == input_mudata.mod["prot"].n_obs
    ), "No prot obs should have been filtered"
    assert (
        mu_out.mod["prot"].n_vars == input_mudata.mod["prot"].n_vars
    ), "No prot vars should have been filtered"
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ], "Feature types of RNA modality should be Gene Expression"
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ], "Feature types of prot modality should be Antibody Capture"


def test_empty_mudata(run_component, random_h5mu_path):
    output_mu = random_h5mu_path()
    empty_mudata_path = random_h5mu_path()
    empty_mudata = mu.MuData(
        {
            modality: ad.AnnData(csr_array((5, 0), dtype=np.int8))
            for modality in ("rna",)
        }
    )

    empty_mudata.write(empty_mudata_path)
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                empty_mudata_path,
                "--output",
                output_mu,
                "--output_compression",
                "gzip",
            ]
        )
    assert re.search(
        "ValueError: Modality rna of input Mudata .* appears to be empty",
        err.value.stdout.decode("utf-8"),
    )


@pytest.mark.xfail(strict=False)
def test_doublet_automatic_threshold_detection_fails(
    run_component, input_with_failed_run, random_h5mu_path
):
    """
    Test if the component fails if doublet score threshold could not automatically be set
    """
    output_mu = random_h5mu_path()

    with pytest.raises(subprocess.CalledProcessError) as e_info:
        run_component(
            [
                "--input",
                input_with_failed_run,
                "--output",
                output_mu,
                "--output_compression",
                "gzip",
                "--num_pca_components",
                "1",
                "--min_gene_variablity_percent",
                "0",
                "--min_cells",
                "1",
                "--min_counts",
                "1",
            ]
        )
    assert re.search(
        r"RuntimeError: Scrublet could not automatically detect the doublet score threshold\. "
        r"--allow_automatic_threshold_detection_fail can be used to ignore this failure and "
        r"set the corresponding output columns to NA\.",
        e_info.value.stdout.decode("utf-8"),
    )

    assert not Path(output_mu).is_file(), "Output file not found"


@pytest.mark.xfail(strict=False)
def test_doublet_automatic_threshold_detection_fails_recovery(
    run_component, input_with_failed_run
):
    """
    Test if the component can recover from scrublet not automatically able to set the doublet score threshold
    and it is not set.
    """
    output_mu = "output-5.h5mu"

    run_component(
        [
            "--input",
            input_with_failed_run,
            "--output",
            output_mu,
            "--output_compression",
            "gzip",
            "--num_pca_components",
            "1",
            "--min_gene_variablity_percent",
            "0",
            "--min_cells",
            "1",
            "--min_counts",
            "1",
            "--allow_automatic_threshold_detection_fail",
        ]
    )
    assert Path(output_mu).is_file(), "Output file not found"

    mu_out = mu.read_h5mu(output_mu)
    assert mu_out.mod["rna"].obs["filter_with_scrublet"].isna().all()


def test_selecting_input_layer(
    run_component, tmp_path, random_h5mu_path, input_mudata, input_mudata_path
):
    output_mu = random_h5mu_path()
    input_data = mu.read_h5mu(input_mudata_path)
    input_data.mod["rna"].layers["test_layer"] = input_data.mod["rna"].X
    input_data.mod["rna"].X = None

    temp_input = tmp_path / "temp.h5mu"
    input_data.write(temp_input)

    run_component(
        [
            "--input",
            temp_input,
            "--output",
            output_mu,
            "--modality",
            "rna",
            "--min_counts",
            "10",
            "--num_pca_components",
            "10",
            "--layer",
            "test_layer",
            "--do_subset",
        ]
    )
    assert Path(output_mu).is_file(), "Output file not found"

    mu_out = mu.read_h5mu(output_mu)
    new_obs = mu_out.mod["rna"].n_obs
    new_vars = mu_out.mod["rna"].n_vars
    assert (
        new_obs < input_mudata.mod["rna"].n_obs
    ), "Some cells should have been filtered"
    assert (
        new_vars == input_mudata.mod["rna"].n_vars
    ), "No genes should have been filtered"
    assert (
        mu_out.mod["prot"].n_obs == input_mudata.mod["prot"].n_obs
    ), "No prot obs should have been filtered"
    assert (
        mu_out.mod["prot"].n_vars == input_mudata.mod["prot"].n_vars
    ), "No prot vars should have been filtered"
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ], "Feature types of RNA modality should be Gene Expression"
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ], "Feature types of prot modality should be Antibody Capture"


if __name__ == "__main__":
    exit(pytest.main([__file__]))
