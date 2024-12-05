from pathlib import Path
import re
import subprocess
import pytest

import mudata as mu
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

## VIASH START
meta = {
    "name": "foo",
    "resources_dir": "resources_test/",
    "executable": "target/executable/filter/filter_with_scrublet/filter_with_scrublet",
}
# def run_component(args_as_list):
#     try:
#         subprocess_args = [meta['executable']] + args_as_list
#         print(" ".join(subprocess_args), flush=True)
#         subprocess.check_output(subprocess_args, stderr=subprocess.STDOUT)
#     except subprocess.CalledProcessError as e:
#         print(e.stdout.decode("utf-8"), flush=True)
#         raise e
## VIASH END

# read input file
input_path = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"
input_mu = mu.read_h5mu(input_path)
orig_obs = input_mu.mod["rna"].n_obs
orig_vars = input_mu.mod["rna"].n_vars
orig_prot_obs = input_mu.mod["prot"].n_obs
orig_prot_vars = input_mu.mod["prot"].n_vars


def test_filter_a_little_bit(run_component):
    output_mu = "output-1.h5mu"

    run_component(
        [
            "--input",
            input_path,
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
    assert new_obs == orig_obs, "No RNA obs should have been filtered"
    assert new_vars == orig_vars, "No RNA vars should have been filtered"
    assert (
        mu_out.mod["prot"].n_obs == orig_prot_obs
    ), "No prot obs should have been filtered"
    assert (
        mu_out.mod["prot"].n_vars == orig_prot_vars
    ), "No prot vars should have been filtered"
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ], "Feature types of RNA modality should be Gene Expression"
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ], "Feature types of prot modality should be Antibody Capture"


def test_filtering_a_lot(run_component):
    output_mu = "output-2.h5mu"

    run_component(
        [
            "--input",
            input_path,
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
    assert new_obs < orig_obs, "Some cells should have been filtered"
    assert new_vars == orig_vars, "No genes should have been filtered"
    assert mu_out.mod["prot"].n_obs == orig_obs, "No prot obs should have been filtered"
    assert (
        mu_out.mod["prot"].n_vars == orig_prot_vars
    ), "No prot vars should have been filtered"
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ], "Feature types of RNA modality should be Gene Expression"
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ], "Feature types of prot modality should be Antibody Capture"


@pytest.fixture(scope="module")
def input_with_failed_run():
    new_mudata_path = "pbmc-perturbed.h5mu"

    mudata_in = mu.read_h5mu(input_path)

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


@pytest.mark.xfail(strict=False)
def test_doublet_automatic_threshold_detection_fails(
    run_component, input_with_failed_run
):
    """
    Test if the component fails if doublet score threshold could not automatically be set
    """
    output_mu = "output-4.h5mu"

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


def test_selecting_input_layer(run_component, tmp_path):
    output_mu = "output-2.h5mu"
    input_data = mu.read_h5mu(input_path)
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
    assert new_obs < orig_obs, "Some cells should have been filtered"
    assert new_vars == orig_vars, "No genes should have been filtered"
    assert mu_out.mod["prot"].n_obs == orig_obs, "No prot obs should have been filtered"
    assert (
        mu_out.mod["prot"].n_vars == orig_prot_vars
    ), "No prot vars should have been filtered"
    assert list(mu_out.mod["rna"].var["feature_types"].cat.categories) == [
        "Gene Expression"
    ], "Feature types of RNA modality should be Gene Expression"
    assert list(mu_out.mod["prot"].var["feature_types"].cat.categories) == [
        "Antibody Capture"
    ], "Feature types of prot modality should be Antibody Capture"


if __name__ == "__main__":
    exit(pytest.main([__file__]))
