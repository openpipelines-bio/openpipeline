import sys
import pytest
import mudata as mu
import numpy as np
from scipy.sparse import issparse

## VIASH START
meta = {"name": "bpcells_regress_out", "resources_dir": "resources_test/"}
## VIASH END


@pytest.fixture
def input_h5mu_path():
    return f"{meta['resources_dir']}/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def output_h5mu_path(tmp_path):
    return tmp_path / "output.h5mu"


def test_regress_out(run_component, input_h5mu_path, output_h5mu_path):
    # execute command
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--output_compression",
        "gzip",
    ]
    run_component(cmd_pars)

    assert output_h5mu_path.is_file(), "No output was created."

    mu_input = mu.read_h5mu(input_h5mu_path)
    mu_output = mu.read_h5mu(output_h5mu_path)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'
    assert "prot" in mu_output.mod, 'Output should contain data.mod["prot"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]
    prot_in = mu_input.mod["prot"]
    prot_out = mu_output.mod["prot"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert prot_in.shape == prot_out.shape, "Should have same shape as before"

    assert np.mean(rna_in.X) != np.mean(rna_out.X), "RNA expression should have changed"
    assert np.mean(prot_in.X) == np.mean(prot_out.X), (
        "Protein expression should remain the same"
    )


def test_no_regress_out_without_obs_keys(
    run_component, input_h5mu_path, output_h5mu_path
):
    # execute command
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
    ]
    run_component(cmd_pars)

    mu_input = mu.read_h5mu(input_h5mu_path)
    mu_output = mu.read_h5mu(output_h5mu_path)

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]

    assert np.mean(rna_in.X) == np.mean(rna_out.X), (
        "RNA expression should remain the same"
    )


def test_regress_out_with_layers(run_component, input_h5mu_path, output_h5mu_path):
    # execute command
    cmd_pars = [
        "--input",
        input_h5mu_path,
        "--output",
        output_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--input_layer",
        "log_normalized",
        "--output_layer",
        "output",
    ]
    run_component(cmd_pars)

    rna_in = mu.read_h5ad(input_h5mu_path, mod="rna")
    rna_out = mu.read_h5ad(output_h5mu_path, mod="rna")

    assert np.mean(rna_in.layers["log_normalized"]) != np.mean(
        rna_out.layers["output"]
    ), "RNA expression should have changed"


def test_regress_out_hvg(run_component, input_h5mu_path, output_h5mu_path, tmp_path):
    base_pars = [
        "--input",
        input_h5mu_path,
        "--obs_keys",
        "total_counts",
        "--input_layer",
        "log_normalized",
        "--output_layer",
        "output",
    ]
    run_component(
        base_pars + ["--output", output_h5mu_path, "--var_input", "filter_with_hvg"]
    )
    all_genes_path = tmp_path / "all_genes.h5mu"
    run_component(base_pars + ["--output", all_genes_path])

    rna_in = mu.read_h5mu(input_h5mu_path).mod["rna"]
    rna_out = mu.read_h5mu(output_h5mu_path).mod["rna"]
    rna_all_genes = mu.read_h5mu(all_genes_path).mod["rna"]
    hvg = rna_in.var["filter_with_hvg"].to_numpy()

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"

    output_layer = rna_out.layers["output"]
    assert issparse(output_layer), "Output layer should be a sparse matrix"
    output_matrix = output_layer.toarray()

    assert not np.any(output_matrix[:, ~hvg]), "Non-selected genes should be set to 0"
    assert output_layer.nnz <= rna_in.n_obs * hvg.sum(), (
        "Only values of selected genes should be stored"
    )
    np.testing.assert_allclose(
        output_matrix[:, hvg],
        rna_all_genes.layers["output"].toarray()[:, hvg],
        err_msg="Selected genes should be regressed as when using all genes",
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
